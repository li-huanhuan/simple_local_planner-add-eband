#include <simple_local_planner/pure_planner_ros.h>
#include <sys/time.h>
#include <boost/tokenizer.hpp>
#include <Eigen/Core>
#include <cmath>
#include <ros/console.h>
#include <pluginlib/class_list_macros.h>
#include <nav_msgs/Path.h>
#include <std_msgs/Bool.h>

//register this planner as a BaseLocalPlanner plugin
PLUGINLIB_EXPORT_CLASS(simple_local_planner::PurePlannerROS, nav_core::BaseLocalPlanner)

namespace simple_local_planner
{

  void PurePlannerROS::reconfigureCB(SimpleLocalPlannerConfig &config, uint32_t level)
  {
      tc_->reconfigure(config);
      eband_->reconfigure(config);
  }

  PurePlannerROS::PurePlannerROS() :
    world_model_(nullptr), tc_(nullptr), dp_(nullptr),costmap_ros_(nullptr), tf_(nullptr), initialized_(false), odom_helper_("odom") {}

  PurePlannerROS::PurePlannerROS(std::string name, tf2_ros::Buffer *tf, costmap_2d::Costmap2DROS *costmap_ros):
    world_model_(nullptr), tc_(nullptr), dp_(nullptr),costmap_ros_(nullptr), tf_(nullptr), initialized_(false), odom_helper_("odom")
  {
    initialize(name, tf, costmap_ros);
  }

  PurePlannerROS::~PurePlannerROS()
  {
    if(dsrv_ != nullptr)
    {
      delete dsrv_;
      dsrv_ = nullptr;
    }

    if(tc_ != nullptr)
    {
      delete tc_;
      tc_ = nullptr;
    }

    if(world_model_ != nullptr)
    {
      delete world_model_;
      world_model_ = nullptr;
    }

    if(dp_ != nullptr)
    {
      delete dp_;
      dp_= nullptr;
    }
  }

  void PurePlannerROS::initialize(std::string name,tf2_ros::Buffer* tf,costmap_2d::Costmap2DROS* costmap_ros)
  {
    if (!isInitialized())
    {
        ros::NodeHandle private_nh("~/" + name);
        l_plan_pub_ = private_nh.advertise<nav_msgs::Path>("local_plan", 1);
        g_replan_pub_ = private_nh.advertise<std_msgs::Bool>("/move_base/replan", 1);

        tf_ = tf;
        costmap_ros_ = costmap_ros;

        //initialize the copy of the costmap the controller will use
        costmap_ = costmap_ros_->getCostmap();

        global_frame_ = costmap_ros_->getGlobalFrameID();
        robot_base_frame_ = costmap_ros_->getBaseFrameID();

        // Robot Configuration Parameters
        private_nh.param("max_trans_vel", max_vel_x_, 0.6);
        private_nh.param("min_trans_vel", min_vel_x_, 0.1);
        private_nh.param("max_rot_vel", max_vel_th_, 0.5);
        private_nh.param("min_rot_vel", min_vel_th_, 0.1);
        private_nh.param("max_trans_acc", max_trans_acc_, 1.0);
        private_nh.param("max_rot_acc", max_rot_acc_, 1.0);
        private_nh.param("min_in_place_rot_vel", min_in_place_vel_th_, 0.3);

        //Goal tolerance parameters
        private_nh.param("yaw_goal_tolerance", yaw_goal_tolerance_, 0.05);
        private_nh.param("xy_goal_tolerance", xy_goal_tolerance_, 0.10);
        private_nh.param("wp_tolerance", wp_tolerance_, 0.5);

        private_nh.param("sim_time", sim_time_, 1.0);
        private_nh.param("sim_granularity", sim_granularity_, 0.025);
        private_nh.param("angular_sim_granularity", angular_sim_granularity_, sim_granularity_);
        private_nh.param("controller_freq", controller_freq_, 15.0);

        bool is_use_dwa;
        private_nh.param("use_dwa", is_use_dwa, false);
        world_model_ = new CostmapModel(*costmap_);
        footprint_spec_ = costmap_ros_->getRobotFootprint();
        tc_ = new PurePlanner(*world_model_, *costmap_, footprint_spec_, controller_freq_,
                              max_vel_x_, min_vel_x_, max_vel_th_, min_vel_th_,
                              min_in_place_vel_th_, max_trans_acc_, max_rot_acc_,
                              yaw_goal_tolerance_, xy_goal_tolerance_, wp_tolerance_,
                              sim_time_, sim_granularity_, angular_sim_granularity_, is_use_dwa);

        // dijkstra planner
        dp_ = new DijkstraPlanner(name,costmap_,global_frame_);
        // eband planner
        eband_ = boost::shared_ptr<EBandPlanner>(new EBandPlanner(name, costmap_ros_));
        // create object for visualization
        eband_visual_ = boost::shared_ptr<EBandVisualization>(new EBandVisualization);
        // pass visualization object to elastic band
        eband_->setVisualization(eband_visual_);
        // initialize visualization - set node handle and pointer to costmap
        eband_visual_->initialize(private_nh, costmap_ros_);

        initialized_ = true;
        //this will load the values of cfg params overwritting the read ones from the yaml file.
        dsrv_ = new dynamic_reconfigure::Server<SimpleLocalPlannerConfig>(private_nh);
        dynamic_reconfigure::Server<SimpleLocalPlannerConfig>::CallbackType cb = boost::bind(&PurePlannerROS::reconfigureCB, this, _1, _2);
        dsrv_->setCallback(cb);
    }
    else
    {
        ROS_WARN("This planner has already been initialized, doing nothing");
    }
  }

  void PurePlannerROS::replanGlobalpath()
  {
    std_msgs::Bool data;
    data.data = true;
    this->g_replan_pub_.publish(data);
  }

  bool PurePlannerROS::checkLocalPath(const std::vector<geometry_msgs::PoseStamped> &path)
  {
    unsigned int mx = 0,my = 0;
    size_t path_size = path.size();
    for (size_t i=0;i<path_size;i++)
    {
      this->costmap_->worldToMap(path[i].pose.position.x,path[i].pose.position.y,mx,my);
      unsigned char cost = this->costmap_->getCost(mx,my);
      if(cost != costmap_2d::FREE_SPACE)
      {
        return false;
      }
    }
    return true;
  }

  bool PurePlannerROS::setPlan(const std::vector<geometry_msgs::PoseStamped>& orig_global_plan)
  {
    std::cout << "received new global plan point size: " << orig_global_plan.size() << std::endl;
    if (! isInitialized())
    {
      ROS_ERROR("This planner has not been initialized, please call initialize() before using this planner");
      return false;
    }

    if(orig_global_plan.empty())
    {
      ROS_ERROR("local planner set plan zero length");
      return false;
    }

    if(!global_plan_.empty())
    {
      double goal_goal_dist =  this->computeDistance(global_plan_.back().pose.position.x,
                                                     global_plan_.back().pose.position.y,
                                                     orig_global_plan.back().pose.position.x,
                                                     orig_global_plan.back().pose.position.y); //新旧目标点的距离
      double goal_goal_aw_diff = this->tc_->normalizeAngle( (tf::getYaw(global_plan_.back().pose.orientation) -
                                                             tf::getYaw(orig_global_plan.back().pose.orientation)),
                                                            -M_PI,
                                                            M_PI); //新旧目标点角度差

      if(goal_goal_dist > 0.1 || goal_goal_aw_diff > 0.1)
      {
        local_plan_.clear();
      }
    }

    //reset the global plan
    global_plan_.clear();
    global_plan_ = orig_global_plan;

    { //add eband param
        std::vector<int> start_end_counts (2, static_cast<int>( global_plan_.size() ) ); // counts from the end() of the plan
        if(!simple_local_planner::transformGlobalPlan(*tf_,
                                                      global_plan_,
                                                      *costmap_ros_,
                                                      costmap_ros_->getGlobalFrameID(),
                                                      local_plan_,
                                                      start_end_counts))
        {
          // if plan could not be tranformed abort control and local planning
          ROS_WARN("Could not transform the global plan to the frame of the controller");
          return false;
        }

        // also check if there really is a plan
        if(local_plan_.empty())
        {
          // if global plan passed in is empty... we won't do anything
          ROS_WARN("Transformed plan is empty. Aborting local planner!");
          return false;
        }

        // set plan - as this is fresh from the global planner robot pose should be identical to start frame
        // 设定计划-因为这是从全局计划程序中重新获得的，所以机器人的姿势应该与开始帧相同。
        if(!eband_->setPlan(local_plan_))
        {
          // We've had some difficulty where the global planner keeps returning a valid path that runs through an obstacle
          // 全球计划人员不断返回一条穿过障碍的有效路径时，我们遇到了一些困难。
          // in the local costmap. See issue #5. Here we clear the local costmap and try one more time.
          costmap_ros_->resetLayers();
          if (!eband_->setPlan(local_plan_))
          {
            ROS_ERROR("Setting plan to Elastic Band method failed!");
            return false;
          }
        }

        // plan transformed and set to elastic band successfully - set counters to global variable
        // 计划已转换并成功设置为弹性带-将计数器设置为全局变量。
        plan_start_end_counter_ = start_end_counts; // 存储转换后的路径在全局计划中的起始和结束编号。

        // let eband refine the plan before starting continuous operation (to smooth sampling based plans)
        // 让eband在开始连续运行之前完善计划（以平滑基于采样的计划）。
        eband_->optimizeBand();
    }

    //reset the goal flag
    reached_goal_ = false;
    return true;
  }

  bool PurePlannerROS::transformGlobalPlan(const tf2_ros::Buffer& tf,
                                           const std::vector<geometry_msgs::PoseStamped>& global_plan,
                                           const geometry_msgs::PoseStamped& global_pose,
                                           const costmap_2d::Costmap2D& costmap,
                                           const std::string& global_frame,
                                           std::vector<geometry_msgs::PoseStamped>& transformed_plan)
  {
    transformed_plan.clear();
    if (global_plan.empty())
    {
      ROS_ERROR("Received plan with zero length");
      return false;
    }
    const geometry_msgs::PoseStamped& plan_pose = global_plan[0];
    try
    {
      // get plan_to_global_transform from plan frame to global_frame
      geometry_msgs::TransformStamped plan_to_global_transform = tf.lookupTransform(global_frame,
                                                                                    ros::Time(),
                                                                                    plan_pose.header.frame_id,
                                                                                    plan_pose.header.stamp,
                                                                                    plan_pose.header.frame_id,
                                                                                    ros::Duration(0.5));

      //let's get the pose of the robot in the frame of the plan
      geometry_msgs::PoseStamped robot_pose;
      tf.transform(global_pose, robot_pose, plan_pose.header.frame_id);

      //we'll discard points on the plan that are outside the local costmap
      double dist_threshold = std::max(costmap.getSizeInCellsX() * costmap.getResolution() / 2.0,
                                       costmap.getSizeInCellsY() * costmap.getResolution() / 2.0) * 0.9;

      unsigned int i = 0;
      double sq_dist_threshold = dist_threshold * dist_threshold;
      double sq_dist = 0;

      //we need to loop to a point on the plan that is within a certain distance of the robot
      unsigned int plan_pose_size = static_cast<unsigned int>(global_plan.size());
      while(i < plan_pose_size)
      {
        double x_diff = robot_pose.pose.position.x - global_plan[i].pose.position.x;
        double y_diff = robot_pose.pose.position.y - global_plan[i].pose.position.y;
        sq_dist = x_diff * x_diff + y_diff * y_diff;
        if (sq_dist <= sq_dist_threshold)
        {
          break;
        }
        ++i;
      }

      geometry_msgs::PoseStamped newer_pose;
      //now we'll transform until points are outside of our distance threshold
      while(i < plan_pose_size && sq_dist <= sq_dist_threshold)
      {
        const geometry_msgs::PoseStamped& pose = global_plan[i];
        tf2::doTransform(pose, newer_pose, plan_to_global_transform);
        transformed_plan.push_back(newer_pose);
        double x_diff = robot_pose.pose.position.x - global_plan[i].pose.position.x;
        double y_diff = robot_pose.pose.position.y - global_plan[i].pose.position.y;
        sq_dist = x_diff * x_diff + y_diff * y_diff;
        ++i;
      }
    }
    catch(tf2::LookupException& ex)
    {
      ROS_ERROR("No Transform available Error: %s\n", ex.what());
      return false;
    }
    catch(tf2::ConnectivityException& ex)
    {
      ROS_ERROR("Connectivity Error: %s\n", ex.what());
      return false;
    }
    catch(tf2::ExtrapolationException& ex)
    {
      ROS_ERROR("Extrapolation Error: %s\n", ex.what());
      if (!global_plan.empty())
        ROS_ERROR("Global Frame: %s Plan Frame size %d: %s\n",
                  global_frame.c_str(),
                  static_cast<unsigned int>(global_plan.size()),
                  global_plan[0].header.frame_id.c_str());

      return false;
    }

    return true;
  }

  void PurePlannerROS::prunePlan(const geometry_msgs::PoseStamped& global_pose, std::vector<geometry_msgs::PoseStamped>& plan, std::vector<geometry_msgs::PoseStamped>& global_plan)
  {
    ROS_ASSERT(global_plan.size() >= plan.size());
    std::vector<geometry_msgs::PoseStamped>::iterator it = plan.begin();
    std::vector<geometry_msgs::PoseStamped>::iterator global_it = global_plan.begin();
    while(it != plan.end())
    {
      const geometry_msgs::PoseStamped& w = *it;
      // Fixed error bound of 2 meters for now. Can reduce to a portion of the map size or based on the resolution
      double x_diff = global_pose.pose.position.x - w.pose.position.x;
      double y_diff = global_pose.pose.position.y - w.pose.position.y;
      double distance_sq = x_diff * x_diff + y_diff * y_diff;
      if(distance_sq < 1)
      {
        break;
      }
      it = plan.erase(it);
      global_it = global_plan.erase(global_it);
    }
  }

  void PurePlannerROS::prunePlan(const geometry_msgs::PoseStamped &robot_pose, std::vector<geometry_msgs::PoseStamped> &plan)
  {
    if(plan.empty())
      return;
    std::vector<geometry_msgs::PoseStamped>::iterator it = plan.begin();
    while(it != plan.end())
    {
      const geometry_msgs::PoseStamped& w = *it;
      // Fixed error bound of 2 meters for now. Can reduce to a portion of the map size or based on the resolution
      double x_diff = robot_pose.pose.position.x - w.pose.position.x;
      double y_diff = robot_pose.pose.position.y - w.pose.position.y;
      double distance_sq = x_diff * x_diff + y_diff * y_diff;
      if(distance_sq < 1)
      {
        break;
      }
      it = plan.erase(it);
    }
  }

  bool PurePlannerROS::computeVelocityCommands(geometry_msgs::Twist& cmd_vel)
  {
    if (! isInitialized())
    {
      ROS_ERROR("This planner has not been initialized, please call initialize() before using this planner");
      return false;
    }

    geometry_msgs::PoseStamped robot_pose;
    if (!costmap_ros_->getRobotPose(robot_pose))
    {
      return false;
    }



    std::vector<geometry_msgs::PoseStamped> tmp_plan;
    // convert robot pose to frame in plan and set position in band at which to append
    // 将机器人的位姿转换成路径中的坐标系，并设置附加到松紧带中的位置。
    tmp_plan.assign(1, robot_pose);
    simple_local_planner::AddAtPosition add_frames_at = add_front;

    // set it to elastic band and let eband connect it
    // 将其设置为松紧带，然后让松紧带连接它。
    if(!eband_->addFrames(tmp_plan, add_frames_at))
    {
      ROS_WARN("Could not connect robot pose to existing elastic band.");
      return false;
    }

    // get additional path-frames which are now in moving window
    // 获取其他路径框架，它们现在正在移动窗口中。
    // Checking for new path frames in moving window.
    // 在移动窗口中检查新的路径框架。
    std::vector<int> plan_start_end_counter = plan_start_end_counter_;
    std::vector<geometry_msgs::PoseStamped> append_transformed_plan;
    // transform global plan to the map frame we are working in - careful this also cuts the plan off (reduces it to local window)
    if(!simple_local_planner::transformGlobalPlan(*tf_,
                                                 global_plan_,
                                                 *costmap_ros_,
                                                 costmap_ros_->getGlobalFrameID(),
                                                 local_plan_,
                                                 plan_start_end_counter))
    {
      // if plan could not be tranformed abort control and local planning
      ROS_WARN("Could not transform the global plan to the frame of the controller");
      return false;
    }

    // also check if there really is a plan
    if(local_plan_.empty())
    {
      // if global plan passed in is empty... we won't do anything
      ROS_WARN("Transformed plan is empty. Aborting local planner!");
      return false;
    }

    // identify new frames - if there are any
    append_transformed_plan.clear(); //(0)存储的为开始点距离全局路径最后一个点的点数差，(1)存储的为转换后的路径里最后一个点与全局路径里最后一个路径点的点数差
    // did last transformed plan end futher away from end of complete plan than this transformed plan?
    if(plan_start_end_counter_.at(1) > plan_start_end_counter.at(1)) // 新转换的计划比旧的转换后的计划更靠近完整计划的末端吗 // counting from the back (as start might be pruned)
    {
      // new frames in moving window
      if(plan_start_end_counter_.at(1) > plan_start_end_counter.at(0)) // counting from the back (as start might be pruned)
      {
        // append everything
        // 追加所有内容。
        append_transformed_plan = local_plan_;
      }
      else
      {
        // append only the new portion of the plan
        // 仅追加计划的新部分。
        int discarded_frames = plan_start_end_counter.at(0) - plan_start_end_counter_.at(1);
        ROS_ASSERT(local_plan_.begin() + discarded_frames + 1 >= local_plan_.begin());
        ROS_ASSERT(local_plan_.begin() + discarded_frames + 1 < local_plan_.end());
        append_transformed_plan.assign(local_plan_.begin() + discarded_frames + 1, local_plan_.end());
      }

      // set it to elastic band and let eband connect it
      // 将其设置为松紧带，然后让松紧带连接它。
      // "Adding %d new frames to current band", (int) append_transformed_plan.size()
      // 将％d个新点添加到当前松紧带。
      if(eband_->addFrames(append_transformed_plan, add_back))
      {
        // appended frames succesfully to global plan - set new start-end counts
        // 成功将框架添加到全局计划-设置新的开始-结束计数。
        // ROS_DEBUG("Sucessfully added frames to band");
        plan_start_end_counter_ = plan_start_end_counter;
      }
      else
      {
        ROS_WARN("Failed to add frames to existing band");
        return false;
      }
    }
    else
    {
      // ROS_DEBUG("Nothing to add");
    }

    // update Elastic Band (react on obstacle from costmap, ...)
    // 更新松紧带（对成本图的障碍物做出反应，...）。
    // Calling optimization method for elastic band.
    // 松紧带的调用优化方法。
    std::vector<simple_local_planner::Bubble> current_band;
    if(!eband_->optimizeBand())
    {
      ROS_WARN("Optimization failed - Band invalid - No controls availlable");
      return false;
    }

    // display current band
    if(eband_->getBand(current_band))
      eband_visual_->publishBand("bubbles", current_band);

    eband_->getPlan(local_plan_);

/**********************************
    double distance_local_goal = 0;
    if(!local_plan_.empty())
      distance_local_goal = computeDistance(robot_pose.pose.position.x,
                                            robot_pose.pose.position.y,
                                            local_plan_.back().pose.position.x,
                                            local_plan_.back().pose.position.y);
    if(distance_local_goal < 1.0)
    {
      if(local_plan_.empty())
      {
        //get the global plan in our frame
        if (!this->transformGlobalPlan(*tf_, global_plan_, robot_pose, *costmap_, global_frame_, local_plan_))
        {
          ROS_WARN("Could not transform the global plan to the frame of the controller");
          return false;
        }
      }
      else if( computeDistance(global_plan_.back().pose.position.x,
                               global_plan_.back().pose.position.y,
                               local_plan_.back().pose.position.x,
                               local_plan_.back().pose.position.y) > 0.025)
      {
        //get the global plan in our frame
        if (!this->transformGlobalPlan(*tf_, global_plan_, robot_pose, *costmap_, global_frame_, local_plan_))
        {
          ROS_WARN("Could not transform the global plan to the frame of the controller");
          return false;
        }
      }
    }

    this->prunePlan(robot_pose, local_plan_);

    //if the global plan passed in is empty... we won't do anything
    size_t plan_pose_size = local_plan_.size();
    if(plan_pose_size == 0)
      return false;

    if(replan_timer<this->controller_freq_*10)
      replan_timer++;

    if(!this->checkLocalPath(local_plan_))
    {
      this->dp_->makePlan(robot_pose,local_plan_[plan_pose_size-1],local_plan_);
      if( replan_timer > this->controller_freq_*2)
      {
        this->replanGlobalpath();
        replan_timer = 0;
      }
    }

**************************/

    if(local_plan_.empty())
      return false;

    //获得机器人当前速度
    tf::Stamped<tf::Pose> robot_vel;
    odom_helper_.getRobotVel(robot_vel);

    //将局部路径传给速度计算模块
    tc_->updatePlan(local_plan_);

    //根据当前机器人速度和局部路径计算速度指令
    bool ok = tc_->findBestAction(robot_pose, robot_vel, cmd_vel);

    //发布局部路径
    this->publishPlan(local_plan_, l_plan_pub_);
    if(!ok)
    {
      ROS_DEBUG_NAMED("trajectory_planner_ros","The rollout planner failed to find a valid plan. This means that the footprint of the robot was in collision for all simulated trajectories.");
      return false;
    }
    return true;
  }

  void PurePlannerROS::publishPlan(const std::vector<geometry_msgs::PoseStamped>& path, const ros::Publisher& pub)
  {
    //given an empty path we won't do anything
    if(path.empty())
      return;

    //create a path message
    nav_msgs::Path gui_path;
    gui_path.poses.resize(path.size());
    gui_path.header.frame_id = path[0].header.frame_id;
    gui_path.header.stamp = path[0].header.stamp;

    // Extract the plan in world co-ordinates, we assume the path is all in the same frame
    for(unsigned int i=0; i < path.size(); i++)
    {
      gui_path.poses[i] = path[i];
    }

    pub.publish(gui_path);
  }

  bool PurePlannerROS::isGoalReached()
  {
    if (! isInitialized())
    {
      ROS_ERROR("This planner has not been initialized, please call initialize() before using this planner");
      return false;
    }
    return tc_->isGoalReached();
  }

};
