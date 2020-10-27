#ifndef PURE_PLANNER_ROS_H
#define PURE_PLANNER_ROS_H

#include <ros/ros.h>
#include <costmap_2d/costmap_2d.h>
#include <costmap_2d/costmap_2d_publisher.h>
#include <costmap_2d/costmap_2d_ros.h>
#include <simple_local_planner/world_model.h>
#include <simple_local_planner/costmap_model.h>
#include <simple_local_planner/pure_planner.h>
#include <tf/transform_datatypes.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/Twist.h>
#include <geometry_msgs/Point.h>
#include <tf/transform_listener.h>
#include <boost/thread.hpp>
#include <string>
#include <angles/angles.h>
#include <nav_core/base_local_planner.h>

//Dynamic reconfigure
#include <dynamic_reconfigure/server.h>
#include <simple_local_planner/SimpleLocalPlannerConfig.h>
#include <simple_local_planner/odometry_helper_ros.h>

#include <dijkstra_planner.h>
#include <eband_local_planner.h>

namespace simple_local_planner {
  /**
   * @class TrajectoryPlannerROS
   * @brief A ROS wrapper for the trajectory controller that queries the param server to construct a controller
   */
  class PurePlannerROS: public nav_core::BaseLocalPlanner
  {
    public:
      /**
       * @brief  Default constructor for the ros wrapper
       */
      PurePlannerROS();

      /**
       * @brief  Constructs the ros wrapper
       * @param name The name to give this instance of the trajectory planner
       * @param tf A pointer to a transform listener
       * @param costmap The cost map to use for assigning costs to trajectories
       */
      PurePlannerROS(std::string name,tf2_ros::Buffer* tf,costmap_2d::Costmap2DROS* costmap_ros);

      /**
       * @brief  Constructs the ros wrapper
       * @param name The name to give this instance of the trajectory planner
       * @param tf A pointer to a transform listener
       * @param costmap The cost map to use for assigning costs to trajectories
       */

      void initialize(std::string name, tf2_ros::Buffer* tf,costmap_2d::Costmap2DROS* costmap_ros);

      /**
       * @brief  Destructor for the wrapper
       */
      ~PurePlannerROS();

      void replanGlobalpath();
      
      bool checkLocalPath(const std::vector<geometry_msgs::PoseStamped>& path);

      /**
       * @brief  Given the current position, orientation, and velocity of the robot,
       * compute velocity commands to send to the base
       * @param cmd_vel Will be filled with the velocity command to be passed to the robot base
       * @return True if a valid trajectory was found, false otherwise
       */
      bool computeVelocityCommands(geometry_msgs::Twist& cmd_vel);

      void publishPlan(const std::vector<geometry_msgs::PoseStamped>& path, const ros::Publisher& pub);

      /**
       * @brief  Set the plan that the controller is following
       * @param orig_global_plan The plan to pass to the controller
       * @return True if the plan was updated successfully, false otherwise
       */
      bool setPlan(const std::vector<geometry_msgs::PoseStamped>& orig_global_plan);

      inline double computeDistance(double x1,double y1,double x2,double y2)
      {
        return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
      }

      bool transformGlobalPlan(
            const tf2_ros::Buffer& tf,
            const std::vector<geometry_msgs::PoseStamped>& global_plan,
            const geometry_msgs::PoseStamped& global_pose,
            const costmap_2d::Costmap2D& costmap,
            const std::string& global_frame,
            std::vector<geometry_msgs::PoseStamped>& transformed_plan);

      void prunePlan(const geometry_msgs::PoseStamped& global_pose,
                std::vector<geometry_msgs::PoseStamped>& plan,
                std::vector<geometry_msgs::PoseStamped>& global_plan);

      void prunePlan(const geometry_msgs::PoseStamped& global_pose,
                std::vector<geometry_msgs::PoseStamped>& plan);


      /**
       * @brief  Check if the goal pose has been achieved
       * @return True if achieved, false otherwise
       */
      bool isGoalReached();
      
      inline bool isInitialized() {return initialized_;}

      /** @brief Return the inner TrajectoryPlanner object.  Only valid after initialize(). */
      PurePlanner* getPlanner() const { return tc_; }

    private:
      /**
       * @brief Callback to update the local planner's parameters based on dynamic reconfigure
       */
      void reconfigureCB(SimpleLocalPlannerConfig &config, uint32_t level);

      WorldModel* world_model_; ///< @brief The world model that the controller will use
      PurePlanner* tc_; ///< @brief The trajectory controller
      DijkstraPlanner* dp_; ///< @brief dijkstra local planner


      boost::shared_ptr<EBandPlanner> eband_;
      boost::shared_ptr<EBandVisualization> eband_visual_;
      std::vector<int> plan_start_end_counter_;


      costmap_2d::Costmap2DROS* costmap_ros_; ///< @brief The ROS wrapper for the costmap the controller will use
      costmap_2d::Costmap2D* costmap_; ///< @brief The costmap the controller will use
      tf2_ros::Buffer* tf_;
      std::string global_frame_; ///< @brief The frame in which the controller will run
      std::string robot_base_frame_; ///< @brief Used as the base frame id of the robot
      std::vector<geometry_msgs::PoseStamped> global_plan_;
      boost::recursive_mutex odom_lock_;
      double controller_freq_;
      // Robot Configuration Parameters
      double max_vel_x_, min_vel_x_;
      double max_vel_th_, min_vel_th_;
      double max_trans_acc_;
      double max_rot_acc_;
      double min_in_place_vel_th_;
      //Goal tolerance parameters
      double yaw_goal_tolerance_;
  	  double xy_goal_tolerance_;
  	  double wp_tolerance_;
      double sim_time_;
  	  double sim_granularity_;
  	  double angular_sim_granularity_;
      bool initialized_;
      bool reached_goal_;
      uint32_t replan_timer = 0;
      ros::Publisher g_replan_pub_,l_plan_pub_;
      dynamic_reconfigure::Server<SimpleLocalPlannerConfig> *dsrv_;
      simple_local_planner::OdometryHelperRos odom_helper_;
      std::vector<geometry_msgs::Point> footprint_spec_;
      std::vector<geometry_msgs::PoseStamped> local_plan_;
  };
};
#endif
