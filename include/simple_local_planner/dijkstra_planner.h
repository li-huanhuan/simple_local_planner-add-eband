#ifndef DIJKSTRA_PLANNER_H
#define DIJKSTRA_PLANNER_H

#include <ros/ros.h>
#include <costmap_2d/costmap_2d.h>
#include <costmap_2d/cost_values.h>
#include <simple_local_planner/world_model.h>
#include <simple_local_planner/trajectory.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/Twist.h>
#include <tf/transform_datatypes.h>
#include <tf2/LinearMath/Matrix3x3.h>
#include <tf2/utils.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>
#include <nav_msgs/OccupancyGrid.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <angles/angles.h>

// inserting onto the priority blocks
#define push_cur(n)  { if (n>=0 && n<ns_ && !pending_[n] && getCost(costs, n)<lethal_cost_ && currentEnd_<PRIORITYBUFSIZE){ currentBuffer_[currentEnd_++]=n; pending_[n]=true; }}
#define push_next(n) { if (n>=0 && n<ns_ && !pending_[n] && getCost(costs, n)<lethal_cost_ &&    nextEnd_<PRIORITYBUFSIZE){    nextBuffer_[   nextEnd_++]=n; pending_[n]=true; }}
#define push_over(n) { if (n>=0 && n<ns_ && !pending_[n] && getCost(costs, n)<lethal_cost_ &&    overEnd_<PRIORITYBUFSIZE){    overBuffer_[   overEnd_++]=n; pending_[n]=true; }}


namespace simple_local_planner {

#define PRIORITYBUFSIZE 10000
#define POT_HIGH 1.0e10
#define INVSQRT2 0.707106781

  class QuadraticCalculator
  {
      public:
          QuadraticCalculator(int nx, int ny)
          {
            setSize(nx,ny);
          }
          ~QuadraticCalculator(){}

          float calculatePotential(float *potential, unsigned char cost, int n, float prev_potential = 0);

          /**
           * @brief  Sets or resets the size of the map
           * @param nx The x size of the map
           * @param ny The y size of the map
           */
          void setSize(int nx, int ny)
          {
              nx_ = nx;
              ny_ = ny;
              ns_ = nx * ny;
          } /**< sets or resets the size of the map */

      protected:
          inline int toIndex(int x, int y)
          {
              return x + nx_ * y;
          }

          int nx_, ny_, ns_; /**< size of grid, in pixels */
};

  class DijkstraExpansion
  {
      public:
          DijkstraExpansion(QuadraticCalculator* p_calc, int nx, int ny);
          ~DijkstraExpansion();
          bool calculatePotentials(unsigned char* costs,
                                   double start_x, double start_y,
                                   double end_x, double end_y,
                                   int cycles,
                                   float* potential);

          /**
           * @brief  Sets or resets the size of the map
           * @param nx The x size of the map
           * @param ny The y size of the map
           */
          void setSize(int nx, int ny); /**< sets or resets the size of the map */

          void setNeutralCost(unsigned char neutral_cost)
          {
              neutral_cost_ = neutral_cost;
              priorityIncrement_ = 2 * neutral_cost_;
          }

          void setPreciseStart(bool precise){ precise_ = precise; }

          void setSafetyControl(bool param){this->safety_control_ = param;}
      private:

          /**
           * @brief  Updates the cell at index n
           * @param costs The costmap
           * @param potential The potential array in which we are calculating
           * @param n The index to update
           */
          void updateCell(unsigned char* costs, float* potential, int n); /** updates the cell at index n */

          float getCost(unsigned char* costs, int n)
          {
              float c = costs[n];
              if (c < lethal_cost_ - 1 || (unknown_ && c-255<0.01))
              {
                  c = c * factor_ + neutral_cost_;
                  if (c >= lethal_cost_)
                      c = lethal_cost_ - 1;
                  return c;
              }
              return lethal_cost_;
          }

          /** block priority buffers */
          int *buffer1_, *buffer2_, *buffer3_; /**< storage buffers for priority blocks */
          int *currentBuffer_, *nextBuffer_, *overBuffer_; /**< priority buffer block ptrs */
          int currentEnd_, nextEnd_, overEnd_; /**< end points of arrays */
          bool *pending_; /**< pending_ cells during propagation */
          bool precise_;
          bool safety_control_ = true;

          /** block priority thresholds */
          float threshold_; /**< current threshold */
          float priorityIncrement_; /**< priority threshold increment */


  public:

      void setLethalCost(unsigned char lethal_cost)
      {
          lethal_cost_ = lethal_cost;
      }
      void setFactor(float factor)
      {
          factor_ = factor;
      }
      void setHasUnknown(bool unknown)
      {
          unknown_ = unknown;
      }

      void clearEndpoint(unsigned char* costs, float* potential, int gx, int gy, int s)
      {
          int startCell = toIndex(gx, gy);
          for(int i=-s;i<=s;i++)
          {
            for(int j=-s;j<=s;j++)
            {
                int n = startCell+i+nx_*j;
                if(potential[n]<POT_HIGH)
                    continue;
                float c = costs[n]+neutral_cost_;
                float pot = p_calc_->calculatePotential(potential, c, n);
                potential[n] = pot;
            }
          }
      }

  protected:
      inline int toIndex(int x, int y)
      {
          return x + nx_ * y;
      }

      int nx_, ny_, ns_; /**< size of grid, in pixels */
      bool unknown_;
      unsigned char lethal_cost_, neutral_cost_;
      int cells_visited_;
      float factor_;
      QuadraticCalculator* p_calc_;

  };

  class GradientPath
  {
      public:
          GradientPath(QuadraticCalculator* p_calc);
          ~GradientPath();

          void setSize(int xs, int ys);

          //
          // Path construction
          // Find gradient at array points, interpolate path
          // Use step size of pathStep, usually 0.5 pixel
          //
          // Some sanity checks:
          //  1. Stuck at same index position
          //  2. Doesn't get near goal
          //  3. Surrounded by high potentials
          //
          bool getPath(float* potential, double start_x, double start_y, double end_x, double end_y, std::vector<std::pair<float, float> >& path);
      private:
          inline int getNearestPoint(int stc, float dx, float dy)
          {
              int pt = stc + static_cast<int>( round(dx) ) + static_cast<int>( (xs_ * round(dy)) );
              return std::max(0, std::min(xs_ * ys_ - 1, pt));
          }
          float gradCell(float* potential, int n);

          float *gradx_, *grady_; /**< gradient arrays, size of potential array */

          float pathStep_; /**< step size for following gradient */


  public:

      inline int getIndex(int x, int y)
      {
          return x + y * xs_;
      }
      void setLethalCost(unsigned char lethal_cost)
      {
          lethal_cost_ = lethal_cost;
      }
  protected:
      int xs_, ys_;
      unsigned char lethal_cost_;
      QuadraticCalculator* p_calc_;
  };

  enum OrientationMode { NONE, FORWARD, INTERPOLATE, FORWARDTHENINTERPOLATE, BACKWARD, LEFTWARD, RIGHTWARD };

  class OrientationFilter
  {
      public:
          OrientationFilter() : omode_(NONE) {}

          void processPath(const geometry_msgs::PoseStamped& start,
                                   std::vector<geometry_msgs::PoseStamped>& path);

          void setAngleBasedOnPositionDerivative(std::vector<geometry_msgs::PoseStamped>& path, int index);
          void interpolate(std::vector<geometry_msgs::PoseStamped>& path,
                           int start_index, int end_index);

          void setMode(OrientationMode new_mode){ omode_ = new_mode; }

          void setWindowSize(int window_size){ window_size_ = window_size; }
      protected:
          OrientationMode omode_;
          int window_size_;
  };

  class DijkstraPlanner
  {
  public:
    DijkstraPlanner();
    DijkstraPlanner(std::string name, costmap_2d::Costmap2D* costmap, std::string frame_id);
    ~DijkstraPlanner();

    void initialize(std::string name, costmap_2d::Costmap2D* costmap, std::string frame_id);

    bool makePlan(const geometry_msgs::PoseStamped start,
                  const geometry_msgs::PoseStamped goal,
                  std::vector<geometry_msgs::PoseStamped>& plan);

    int optimizationPath(std::vector<geometry_msgs::PoseStamped>& plan,double movement_angle_range);

    void optimizationOrientation(std::vector<geometry_msgs::PoseStamped> &plan);

    double inline distance(double x1,double y1,double x2,double y2);

    bool inline isAroundFree(unsigned int mx,unsigned int my);

    void getNearFreePoint(const geometry_msgs::PoseStamped in,
                          geometry_msgs::PoseStamped& out,
                          double tolerance);

    /**
     * @brief  Computes the full navigation function for the map given a point in the world to start from
     * @param world_point The point to use for seeding the navigation function
     * @return True if the navigation function was computed successfully, false otherwise
     */
    bool computePotential(const geometry_msgs::Point& world_point);

    /**
     * @brief Compute a plan to a goal after the potential for a start point has already been computed (Note: You should call computePotential first)
     * @param start_x
     * @param start_y
     * @param end_x
     * @param end_y
     * @param goal The goal pose to create a plan to
     * @param plan The plan... filled by the planner
     * @return True if a valid plan was found, false otherwise
     */
    bool getPlanFromPotential(double start_x, double start_y,
                              double end_x, double end_y,
                              const geometry_msgs::PoseStamped& goal,
                              std::vector<geometry_msgs::PoseStamped>& plan);

    /**
     * @brief Get the potential, or naviagation cost, at a given point in the world (Note: You should call computePotential first)
     * @param world_point The point to get the potential for
     * @return The navigation function's value at that point in the world
     */
    double getPointPotential(const geometry_msgs::Point& world_point);

    /**
     * @brief Check for a valid potential value at a given point in the world (Note: You should call computePotential first)
     * @param world_point The point to get the potential for
     * @return True if the navigation function is valid at that point in the world, false otherwise
     */
    bool validPointPotential(const geometry_msgs::Point& world_point);

    /**
     * @brief Check for a valid potential value at a given point in the world (Note: You should call computePotential first)
     * @param world_point The point to get the potential for
     * @param tolerance The tolerance on searching around the world_point specified
     * @return True if the navigation function is valid at that point in the world, false otherwise
     */
    bool validPointPotential(const geometry_msgs::Point& world_point, double tolerance);

protected:

    /**
     * @brief Store a copy of the current costmap in \a costmap.  Called by makePlan.
     */
    costmap_2d::Costmap2D* costmap_;
    std::string frame_id_;
    bool initialized_, allow_unknown_;

private:
    void mapToWorld(double mx, double my, double& wx, double& wy);
    bool worldToMap(double wx, double wy, double& mx, double& my);
    void clearRobotCell(const geometry_msgs::PoseStamped& global_pose, unsigned int mx, unsigned int my);
    void publishPotential(float* potential);
    void outlineMap(unsigned char* costarr, int nx, int ny, unsigned char value);
    double inline normalizeAngle(double val,double min = -M_PI,double max = M_PI)
    {
      float norm = 0.0;
      if (val >= min)
        norm = min + fmod((val - min), (max-min));
      else
        norm = max - fmod((min - val), (max-min));

      return norm;
    }

    boost::mutex mutex_;

    QuadraticCalculator* p_calc_;
    DijkstraExpansion* planner_;
    GradientPath* path_maker_;
    OrientationFilter* orientation_filter_;
    ros::Publisher potential_pub_;
    int publish_scale_;
    double default_tolerance_;
    unsigned char* cost_array_;
    unsigned int start_x_, start_y_;
    float* potential_array_;
    float convert_offset_;
  };
}

#endif


