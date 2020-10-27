#include <conversions_and_types.h>

namespace simple_local_planner{

  void PoseToPose2D(const geometry_msgs::Pose pose, geometry_msgs::Pose2D& pose2D)
  {
    // use tf-pkg to convert angles
    tf2::Transform pose_tf;

    // convert geometry_msgs::PoseStamped to tf2::Transform
    tf2::convert(pose, pose_tf);

    // now get Euler-Angles from pose_tf
    double useless_pitch, useless_roll, yaw;
    pose_tf.getBasis().getEulerYPR(yaw, useless_pitch, useless_roll);

    // normalize angle
    yaw = angles::normalize_angle(yaw);

    // and set to pose2D
    pose2D.x = pose.position.x;
    pose2D.y = pose.position.y;
    pose2D.theta = yaw;
  }


  void Pose2DToPose(geometry_msgs::Pose& pose, const geometry_msgs::Pose2D pose2D)
  {
    // use tf-pkg to convert angles
    tf2::Quaternion frame_quat;

    // transform angle from euler-angle to quaternion representation
    frame_quat.setRPY(0.0, 0.0, pose2D.theta);

    // set position
    pose.position.x = pose2D.x;
    pose.position.y = pose2D.y;
    pose.position.z = 0.0;

    // set quaternion
    pose.orientation.x = frame_quat.x();
    pose.orientation.y = frame_quat.y();
    pose.orientation.z = frame_quat.z();
    pose.orientation.w = frame_quat.w();
  }


  bool transformGlobalPlan(const tf2_ros::Buffer& tf,
                           const std::vector<geometry_msgs::PoseStamped>& global_plan,
                           costmap_2d::Costmap2DROS& costmap,
                           const std::string& global_frame,
                           std::vector<geometry_msgs::PoseStamped>& transformed_plan,
                           std::vector<int>& start_end_counts)
  {
    if (global_plan.empty())
    {
      ROS_ERROR("Recieved plan with zero length");
      return false;
    }

    const geometry_msgs::PoseStamped& plan_pose = global_plan[0];

    // initiate refernce variables
    transformed_plan.clear();

    try
    {
      geometry_msgs::TransformStamped transform;
      transform = tf.lookupTransform(global_frame, ros::Time(), plan_pose.header.frame_id, plan_pose.header.stamp, plan_pose.header.frame_id);

      //let's get the pose of the robot in the frame of the plan
      geometry_msgs::TransformStamped robot_pose;
      robot_pose.transform.rotation.w = 1; // identity quaternion
      robot_pose.header.frame_id = costmap.getBaseFrameID();
      robot_pose.header.stamp = ros::Time();
      tf.transform(robot_pose, robot_pose, plan_pose.header.frame_id);
      //we'll keep points on the plan that are within the window that we're looking at

      double dist_threshold = std::max(costmap.getCostmap()->getSizeInMetersX() / 2.0,
                                       costmap.getCostmap()->getSizeInMetersY() / 2.0);

      unsigned int i = 0;
      double sq_dist_threshold = dist_threshold * dist_threshold;
      double sq_dist = DBL_MAX;

      // --- start - modification w.r.t. base_local_planner
      // initiate start_end_count
      std::vector<int> start_end_count;
      start_end_count.assign(2, 0);

      // we know only one direction and that is forward! - initiate search with previous start_end_counts
      // this is neccesserry to work with the sampling based planners - path may severall time enter and leave moving window
      //我们只知道一个方向，那是前进的方向！ -使用以前的start_end_counts开始搜索，这对于与基于采样的计划者一起工作是必要的-路径可能需要数次进入和离开移动窗口。
      ROS_ASSERT( (start_end_counts.at(0) > 0) && (start_end_counts.at(0) <= static_cast<int>( global_plan.size() ) ) );
      i = static_cast<unsigned int>(global_plan.size()) - static_cast<unsigned int>(start_end_counts.at(0));
      // --- end - modification w.r.t. base_local_planner

      //we need to loop to a point on the plan that is within a certain distance of the robot
      while(i < static_cast<unsigned int>(global_plan.size()) && sq_dist > sq_dist_threshold)
      {
       double x_diff = robot_pose.transform.translation.x - global_plan[i].pose.position.x;
       double y_diff = robot_pose.transform.translation.y - global_plan[i].pose.position.y;
       sq_dist = x_diff * x_diff + y_diff * y_diff;

       // --- start - modification w.r.t. base_local_planner
       // not yet in reach - get next frame
       if( sq_dist > sq_dist_threshold )
         ++i;
       else
         // set counter for start of transformed intervall - from back as beginning of plan might be prunned
         start_end_count.at(0) = static_cast<int>( static_cast<unsigned int>(global_plan.size()) - i);
       // --- end - modification w.r.t. base_local_planner
      }


       // tf::Stamped<tf::Pose> tf_pose;
      tf2::Stamped<tf2::Transform> tf_pose, tf_transform;
      geometry_msgs::PoseStamped newer_pose;

      //now we'll transform until points are outside of our distance threshold
      while(i < static_cast<unsigned int>(global_plan.size()) && sq_dist < sq_dist_threshold)
      {
       double x_diff = robot_pose.transform.translation.x - global_plan[i].pose.position.x;
       double y_diff = robot_pose.transform.translation.y - global_plan[i].pose.position.y;
       sq_dist = x_diff * x_diff + y_diff * y_diff;

       const geometry_msgs::PoseStamped& pose = global_plan[i];
       // poseStampedMsgToTF(pose, tf_pose);
       // tf_pose.setData(transform * tf_pose);
       // tf_pose.stamp_ = transform.stamp_;
       tf2::convert(pose, tf_pose);
       tf2::convert(transform, tf_transform);
       tf_pose.setData(tf_transform * tf_pose);
       tf_pose.stamp_ = tf_transform.stamp_;
       tf_pose.frame_id_ = global_frame;
       tf2::toMsg(tf_pose, newer_pose);
       transformed_plan.push_back(newer_pose);

       // --- start - modification w.r.t. base_local_planner
       // set counter for end of transformed intervall - from back as beginning of plan might be prunned
       start_end_count.at(1) = static_cast<int>( static_cast<unsigned int>(global_plan.size()) - i );
       // --- end - modification w.r.t. base_local_planner

       ++i;
      }

       // --- start - modification w.r.t. base_local_planner
       // write to reference variable
       start_end_counts = start_end_count;
      // --- end - modification w.r.t. base_local_planner
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

  double getCircumscribedRadius(costmap_2d::Costmap2DROS& costmap)
  {
    std::vector<geometry_msgs::Point> footprint(costmap.getRobotFootprint());
    double max_distance_sqr = 0;
    for (size_t i = 0; i < footprint.size(); ++i)
    {
      geometry_msgs::Point& p = footprint[i];
      double distance_sqr = p.x*p.x + p.y*p.y;
      if (distance_sqr > max_distance_sqr)
      {
        max_distance_sqr = distance_sqr;
      }
    }
    return sqrt(max_distance_sqr);
  }


}
