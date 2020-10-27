#include <dijkstra_planner.h>

namespace simple_local_planner{

  float QuadraticCalculator::calculatePotential(float *potential, unsigned char cost, int n, float prev_potential)
  {
      prev_potential = 0;

      // get neighbors
      float u, d, l, r;
      l = potential[n - 1];
      r = potential[n + 1];
      u = potential[n - nx_];
      d = potential[n + nx_];
      //  ROS_INFO("[Update] c: %f  l: %f  r: %f  u: %f  d: %f\n",
      //     potential[n], l, r, u, d);
      //  ROS_INFO("[Update] cost: %d\n", costs[n]);

      // find lowest, and its lowest neighbor
      float ta, tc;
      if (l < r)
          tc = l;
      else
          tc = r;
      if (u < d)
          ta = u;
      else
          ta = d;

      float hf = cost; // traversability factor
      float dc = tc - ta;        // relative cost between ta,tc
      if (dc < 0)         // tc is lowest
      {
          dc = -dc;
          ta = tc;
      }

      // calculate new potential
      if (dc >= hf)        // if too large, use ta-only update
          return ta + hf;
      else            // two-neighbor interpolation update
      {
          // use quadratic approximation
          // might speed this up through table lookup, but still have to
          //   do the divide
          float d = dc / hf;
          float v = -0.2301 * d * d + 0.5307 * d + 0.7040;
          return ta + hf * v;
      }
  }

  DijkstraExpansion::DijkstraExpansion(QuadraticCalculator *p_calc, int nx, int ny):
    pending_(nullptr),unknown_(true), lethal_cost_(253), neutral_cost_(50), factor_(3.0), p_calc_(p_calc)
  {
    setSize(nx, ny);
    buffer1_ = new int[PRIORITYBUFSIZE];
    buffer2_ = new int[PRIORITYBUFSIZE];
    buffer3_ = new int[PRIORITYBUFSIZE];
    priorityIncrement_ = 2 * neutral_cost_;
  }

  DijkstraExpansion::~DijkstraExpansion()
  {
    delete[] buffer1_;
    delete[] buffer2_;
    delete[] buffer3_;
    if (pending_)
        delete[] pending_;
  }

  void DijkstraExpansion::setSize(int xs, int ys) //Set/Reset map size
  {
      nx_ = xs;
      ny_ = ys;
      ns_ = xs * ys;
      if (pending_)
          delete[] pending_;

      pending_ = new bool[ns_];
      memset(pending_, 0, ns_ * sizeof(bool));
  }

//
// main propagation function
// Dijkstra method, breadth-first
// runs for a specified number of cycles,
// or until it runs out of cells to update,
// or until the Start cell is found (atStart = true)
//

  bool DijkstraExpansion::calculatePotentials(unsigned char* costs,
                                              double start_x, double start_y,
                                              double end_x, double end_y,
                                              int cycles,
                                              float* potential)
  {
      cells_visited_ = 0;
      // priority buffers
      threshold_ = lethal_cost_;
      currentBuffer_ = buffer1_;
      currentEnd_ = 0;
      nextBuffer_ = buffer2_;
      nextEnd_ = 0;
      overBuffer_ = buffer3_;
      overEnd_ = 0;
      memset(pending_, 0, ns_ * sizeof(bool));
      std::fill(potential, potential + ns_, POT_HIGH);

      // set goal
      int k = toIndex(start_x, start_y);

      if(precise_)
      {
          double dx = start_x - static_cast<int>(start_x), dy = start_y - static_cast<int>(start_y);
          dx = floorf(dx * 100 + 0.5) / 100;
          dy = floorf(dy * 100 + 0.5) / 100;
          potential[k] = neutral_cost_ * 2 * dx * dy;
          potential[k+1] = neutral_cost_ * 2 * (1-dx)*dy;
          potential[k+nx_] = neutral_cost_*2*dx*(1-dy);
          potential[k+nx_+1] = neutral_cost_*2*(1-dx)*(1-dy);//*/

          push_cur(k+2);
          push_cur(k-1);
          push_cur(k+nx_-1);
          push_cur(k+nx_+2);

          push_cur(k-nx_);
          push_cur(k-nx_+1);
          push_cur(k+nx_*2);
          push_cur(k+nx_*2+1);
      }
      else
      {
          potential[k] = 0;
          push_cur(k+1);
          push_cur(k-1);
          push_cur(k-nx_);
          push_cur(k+nx_);
      }

      int nwv = 0;            // max priority block size
      int nc = 0;            // number of cells put into priority blocks
      int cycle = 0;        // which cycle we're on

      // set up start cell
      int startCell = toIndex(end_x, end_y);

      for (; cycle < cycles; cycle++) // go for this many cycles, unless interrupted
      {
          if (currentEnd_ == 0 && nextEnd_ == 0) // priority blocks empty
              return false;
          // stats
          nc += currentEnd_;
          if (currentEnd_ > nwv)
              nwv = currentEnd_;

          // reset pending_ flags on current priority buffer
          int *pb = currentBuffer_;
          int i = currentEnd_;
          while (i-- > 0)
              pending_[*(pb++)] = false;

          // process current priority buffer
          pb = currentBuffer_;
          i = currentEnd_;
          while (i-- > 0)
              updateCell(costs, potential, *pb++);

          // swap priority blocks currentBuffer_ <=> nextBuffer_
          currentEnd_ = nextEnd_;
          nextEnd_ = 0;
          pb = currentBuffer_;        // swap buffers
          currentBuffer_ = nextBuffer_;
          nextBuffer_ = pb;

          // see if we're done with this priority level
          if (currentEnd_ == 0)
          {
              threshold_ += priorityIncrement_;    // increment priority threshold
              currentEnd_ = overEnd_;    // set current to overflow block
              overEnd_ = 0;
              pb = currentBuffer_;        // swap buffers
              currentBuffer_ = overBuffer_;
              overBuffer_ = pb;
          }

          // check if we've hit the Start cell
          if (potential[startCell] < POT_HIGH)
              break;
      }
      if (cycle < cycles)
          return true; // finished up here
      else
          return false;
  }

//
// Critical function: calculate updated potential value of a cell,
// given its neighbors' values
// Planar-wave update calculation from two lowest neighbors in a 4-grid
// Quadratic approximation to the interpolated value
// No checking of bounds here, this function should be fast
//

  inline void DijkstraExpansion::updateCell(unsigned char* costs, float* potential, int n)
  {
      cells_visited_++;

      // do planar wave update
      float c = getCost(costs, n);
      if (c >= lethal_cost_)    // don't propagate into obstacles
          return;

      if(this->safety_control_)
      {
        float c_near[8] = {0};
        int np = n - nx_;
        int nn = n + nx_;
        c_near[0] = getCost(costs, np-1);
        c_near[1] = getCost(costs, np);
        c_near[2] = getCost(costs, np+1);

        c_near[3] = getCost(costs, n-1);
        c_near[4] = getCost(costs, n+1);

        c_near[5] = getCost(costs, nn-1);
        c_near[6] = getCost(costs, nn);
        c_near[7] = getCost(costs, nn+1);

        for(unsigned int i=0;i<8;i++)
        {
          if(c_near[i] > 50)
          {
            return;
          }
        }
      }

      float pot = p_calc_->calculatePotential(potential, c, n);

      // now add affected neighbors to priority blocks
      if (pot < potential[n])
      {
          float le = INVSQRT2 * static_cast<float>( getCost(costs, n - 1) );
          float re = INVSQRT2 * static_cast<float>( getCost(costs, n + 1) );
          float ue = INVSQRT2 * static_cast<float>( getCost(costs, n - nx_) );
          float de = INVSQRT2 * static_cast<float>( getCost(costs, n + nx_) );
          potential[n] = pot;
          if (pot < threshold_)    // low-cost buffer block
          {
              if (potential[n - 1] > pot + le)
                  push_next(n-1);
              if (potential[n + 1] > pot + re)
                  push_next(n+1);
              if (potential[n - nx_] > pot + ue)
                  push_next(n-nx_);
              if (potential[n + nx_] > pot + de)
                  push_next(n+nx_);
          }
          else            // overflow block
          {
              if (potential[n - 1] > pot + le)
                  push_over(n-1);
              if (potential[n + 1] > pot + re)
                  push_over(n+1);
              if (potential[n - nx_] > pot + ue)
                  push_over(n-nx_);
              if (potential[n + nx_] > pot + de)
                  push_over(n+nx_);
          }
      }
  }

  GradientPath::GradientPath(QuadraticCalculator *p_calc) : pathStep_(0.5)
  {
    p_calc_ = p_calc;
    gradx_ = nullptr;
    grady_ = nullptr;
  }

  GradientPath::~GradientPath()
  {
      if (gradx_)
          delete[] gradx_;
      if (grady_)
          delete[] grady_;
  }

  void GradientPath::setSize(int xs, int ys)
  {
      xs_ = xs;
      ys_ = ys;

      if (gradx_)
          delete[] gradx_;
      if (grady_)
          delete[] grady_;
      gradx_ = new float[xs * ys];
      grady_ = new float[xs * ys];
  }

  bool GradientPath::getPath(float* potential, double start_x, double start_y, double goal_x, double goal_y, std::vector<std::pair<float, float> >& path)
  {
      std::pair<float, float> current;
      int stc = getIndex(static_cast<int>(goal_x), static_cast<int>(goal_y));

      // set up offset
      float dx = 0.0; //goal_x - (int)goal_x;
      float dy = 0.0; //goal_y - (int)goal_y;
      int ns = xs_ * ys_;
      memset(gradx_, 0, static_cast<unsigned long>(ns) * sizeof(float));
      memset(grady_, 0, static_cast<unsigned long>(ns) * sizeof(float));

      int c = 0;
      while (c++<ns*4)
      {
          float nx = stc % xs_ + dx, ny = stc / xs_ + dy;
          // check if near goal
          if (fabs(nx - start_x) < 1.5 && fabs(ny - start_y) < 1.5)
          {
              current.first = start_x;
              current.second = start_y;
              path.push_back(current);
              return true;
          }

          if (stc < xs_ || stc > xs_ * ys_ - xs_) // would be out of bounds
          {
              ROS_INFO("would be out of bounds");
              return false;
          }

          current.first = nx;
          current.second = ny;

          path.push_back(current);

          bool oscillation_detected = false;
          int npath = path.size();
          if (npath > 2 && fabs( path[npath - 1].first - path[npath - 3].first ) < 0.0001
                        && fabs( path[npath - 1].second - path[npath - 3].second) < 0.0001 )
          {
              oscillation_detected = true;
          }

          int stcnx = stc + xs_;
          int stcpx = stc - xs_;

          // check for potentials at eight positions near cell
          if (potential[stc] >= POT_HIGH || potential[stc + 1] >= POT_HIGH || potential[stc - 1] >= POT_HIGH
              || potential[stcnx] >= POT_HIGH || potential[stcnx + 1] >= POT_HIGH || potential[stcnx - 1] >= POT_HIGH
              || potential[stcpx] >= POT_HIGH || potential[stcpx + 1] >= POT_HIGH || potential[stcpx - 1] >= POT_HIGH
              || oscillation_detected)
          { //产生了振荡或者8个方向有其中一个为障碍物
              //ROS_DEBUG("[Path] Pot fn boundary, following grid (%0.1f/%d)", potential[stc], (int) path.size());
              // check eight neighbors to find the lowest
              int minc = stc;
              int minp = potential[stc];
              int st = stcpx - 1;
              if (potential[st] < minp)
              {
                  minp = potential[st];
                  minc = st;
              }
              st++;
              if (potential[st] < minp)
              {
                  minp = potential[st];
                  minc = st;
              }
              st++;
              if (potential[st] < minp)
              {
                  minp = potential[st];
                  minc = st;
              }
              st = stc - 1;
              if (potential[st] < minp)
              {
                  minp = potential[st];
                  minc = st;
              }
              st = stc + 1;
              if (potential[st] < minp)
              {
                  minp = potential[st];
                  minc = st;
              }
              st = stcnx - 1;
              if (potential[st] < minp)
              {
                  minp = potential[st];
                  minc = st;
              }
              st++;
              if (potential[st] < minp)
              {
                  minp = potential[st];
                  minc = st;
              }
              st++;
              if (potential[st] < minp)
              {
                  minp = potential[st];
                  minc = st;
              }
              stc = minc;
              dx = 0;
              dy = 0;

              if (potential[stc] >= POT_HIGH)
              {
                  return false;
              }
          }
          else // have a good gradient here 八个方向都没有障碍物
          {
              // get grad at four positions near cell
              gradCell(potential, stc);
              gradCell(potential, stc + 1);
              gradCell(potential, stcnx);
              gradCell(potential, stcnx + 1);

              // get interpolated gradient
              float x1 = (1.0 - dx) * gradx_[stc] + dx * gradx_[stc + 1];
              float x2 = (1.0 - dx) * gradx_[stcnx] + dx * gradx_[stcnx + 1];
              float x = (1.0 - dy) * x1 + dy * x2; // interpolated x
              float y1 = (1.0 - dx) * grady_[stc] + dx * grady_[stc + 1];
              float y2 = (1.0 - dx) * grady_[stcnx] + dx * grady_[stcnx + 1];
              float y = (1.0 - dy) * y1 + dy * y2; // interpolated y

              // check for zero gradient, failed
              if (x == 0.0 && y == 0.0)
              {
                  return false;
              }

              // move in the right direction
              float ss = pathStep_ / hypot(x, y);
              dx += x * ss;
              dy += y * ss;

              // check for overflow
              if (dx > 1.0)
              {
                  stc++;
                  dx -= 1.0;
              }
              if (dx < -1.0)
              {
                  stc--;
                  dx += 1.0;
              }
              if (dy > 1.0)
              {
                  stc += xs_;
                  dy -= 1.0;
              }
              if (dy < -1.0)
              {
                  stc -= xs_;
                  dy += 1.0;
              }
          }
      }
      return false;
  }

// gradient calculations
// calculate gradient at a cell
// positive value are to the right and down
  float GradientPath::gradCell(float* potential, int n)
  {
      if (gradx_[n] + grady_[n] > 0.0)    // check this cell
          return 1.0;

      if (n < xs_ || n > xs_ * ys_ - xs_)    // would be out of bounds
          return 0.0;
      float cv = potential[n];
      float dx = 0.0;
      float dy = 0.0;

      // check for in an obstacle
      if (cv >= POT_HIGH) //如果当前点是障碍物
      {
          if (potential[n - 1] < POT_HIGH) //X方向的前一个点不是障碍物
              dx = -lethal_cost_; //梯度增量很小
          else if (potential[n + 1] < POT_HIGH) //X方向的下一个点不是障碍物
              dx = lethal_cost_; //梯度增量很大

          if (potential[n - xs_] < POT_HIGH) //y方向的上一个点不是障碍物
              dy = -lethal_cost_; //梯度增量很小
          else if (potential[n + xs_] < POT_HIGH) //y方向的下一个点不是障碍物
              dy = lethal_cost_;  //梯度增量很大
      }
      else // not in an obstacle
      {
          // dx calc, average to sides
        //X方向的梯度增量为（X方向上前一个点的势减去当前点的势）+（X方向上当前点的势减去下一个点的势）
          if (potential[n - 1] < POT_HIGH) //X方向的前一个点不是障碍物
              dx += potential[n - 1] - cv; //梯度加上前一个点的势减去当前点的势
          if (potential[n + 1] < POT_HIGH)
              dx += cv - potential[n + 1];

          // dy calc, average to sides
          //Y方向的梯度增量为（Y方向上前一个点的势减去当前点的势）+（Y方向上当前点的势减去下一个点的势）
          if (potential[n - xs_] < POT_HIGH)
              dy += potential[n - xs_] - cv;
          if (potential[n + xs_] < POT_HIGH)
              dy += cv - potential[n + xs_];
      }
      // normalize
      float norm = hypot(dx, dy); //两个方向上取模
      if (norm > 0)
      {
          norm = 1.0 / norm;
          gradx_[n] = norm * dx;
          grady_[n] = norm * dy;
      }
      return norm;
  }

  void set_angle(geometry_msgs::PoseStamped* pose, double angle)
  {
    tf2::Quaternion q;
    q.setRPY(0, 0, angle);
    tf2::convert(q, pose->pose.orientation);
  }

  void OrientationFilter::processPath(const geometry_msgs::PoseStamped& start,
                                      std::vector<geometry_msgs::PoseStamped>& path)
  {
      int n = path.size();
      if (n == 0) return;
      switch(omode_)
      {
          case NONE:
              break;
          case FORWARD:
              for(int i=0;i<n-1;i++)
              {
                  setAngleBasedOnPositionDerivative(path, i);
              }
              break;
          case BACKWARD:
              for(int i=0;i<n-1;i++)
              {
                  setAngleBasedOnPositionDerivative(path, i);
                  set_angle(&path[i], angles::normalize_angle(tf2::getYaw(path[i].pose.orientation) + M_PI));
              }
              break;
          case LEFTWARD:
              for(int i=0;i<n-1;i++)
              {
                  setAngleBasedOnPositionDerivative(path, i);
                  set_angle(&path[i], angles::normalize_angle(tf2::getYaw(path[i].pose.orientation) - M_PI_2));
              }
              break;
          case RIGHTWARD:
              for(int i=0;i<n-1;i++)
              {
                  setAngleBasedOnPositionDerivative(path, i);
                  set_angle(&path[i], angles::normalize_angle(tf2::getYaw(path[i].pose.orientation) + M_PI_2));
              }
              break;
          case INTERPOLATE:
              path[0].pose.orientation = start.pose.orientation;
              interpolate(path, 0, n-1);
              break;
          case FORWARDTHENINTERPOLATE:
              for(int i=0;i<n-1;i++)
              {
                  setAngleBasedOnPositionDerivative(path, i);
              }
              int i=n-3;
              const double last = tf2::getYaw(path[i].pose.orientation);
              while( i>0 )
              {
                  const double new_angle = tf2::getYaw(path[i-1].pose.orientation);
                  double diff = fabs(angles::shortest_angular_distance(new_angle, last));
                  if( diff>0.35)
                      break;
                  else
                      i--;
              }
              path[0].pose.orientation = start.pose.orientation;
              interpolate(path, i, n-1);
              break;
      }
  }

  void OrientationFilter::setAngleBasedOnPositionDerivative(std::vector<geometry_msgs::PoseStamped>& path, int index)
  {
    int index0 = std::max(0, index - window_size_);
    int index1 = std::min( static_cast<int>(path.size() - 1), index + window_size_ );

    double x0 = path[index0].pose.position.x,
           y0 = path[index0].pose.position.y,
           x1 = path[index1].pose.position.x,
           y1 = path[index1].pose.position.y;

    double angle = atan2(y1-y0,x1-x0);
    set_angle(&path[index], angle);
  }

  void OrientationFilter::interpolate(std::vector<geometry_msgs::PoseStamped>& path,
                                      int start_index, int end_index)
  {
      const double start_yaw = tf2::getYaw(path[start_index].pose.orientation),
                   end_yaw   = tf2::getYaw(path[end_index  ].pose.orientation);
      double diff = angles::shortest_angular_distance(start_yaw, end_yaw);
      double increment = diff/(end_index-start_index);
      for(int i=start_index; i<=end_index; i++){
          double angle = start_yaw + increment * i;
          set_angle(&path[i], angle);
      }
  }

  DijkstraPlanner::DijkstraPlanner() :
    costmap_(nullptr), initialized_(false), allow_unknown_(true)
  {
    ;
  }

  DijkstraPlanner::DijkstraPlanner(std::string name, costmap_2d::Costmap2D *costmap, std::string frame_id) :
    costmap_(nullptr), initialized_(false), allow_unknown_(true)
  {
    initialize(name, costmap, frame_id);
  }

  DijkstraPlanner::~DijkstraPlanner()
  {
    if (p_calc_)
        delete p_calc_;
    if (planner_)
        delete planner_;
    if (path_maker_)
        delete path_maker_;
  }

  void DijkstraPlanner::initialize(std::string name, costmap_2d::Costmap2D *costmap, std::string frame_id)
  {
    if (!initialized_)
    {
        ros::NodeHandle private_nh("~/" + name);
        potential_pub_ = private_nh.advertise<nav_msgs::OccupancyGrid>("potential",1);
        costmap_ = costmap;
        frame_id_ = frame_id;
        unsigned int cx = costmap->getSizeInCellsX(), cy = costmap->getSizeInCellsY();
        convert_offset_ = 0.5;

        p_calc_ = new QuadraticCalculator(cx, cy);
        planner_ = new DijkstraExpansion(p_calc_, cx, cy);
        planner_->setPreciseStart(true);
        path_maker_ = new GradientPath(p_calc_);
        orientation_filter_ = new OrientationFilter();

        private_nh.param("allow_unknown", allow_unknown_, false);
        planner_->setHasUnknown(allow_unknown_);
        private_nh.param("default_tolerance", default_tolerance_, 0.0);
        private_nh.param("publish_scale", publish_scale_, 100);

        initialized_ = true;
    }
    else
    {
      ROS_WARN("This planner has already been initialized, you can't call it twice, doing nothing");
    }
  }

  void DijkstraPlanner::outlineMap(unsigned char* costarr, int nx, int ny, unsigned char value)
  {
      unsigned char* pc = costarr;
      for (int i = 0; i < nx; i++)
          *pc++ = value;
      pc = costarr + (ny - 1) * nx;
      for (int i = 0; i < nx; i++)
          *pc++ = value;
      pc = costarr;
      for (int i = 0; i < ny; i++, pc += nx)
          *pc = value;
      pc = costarr + nx - 1;
      for (int i = 0; i < ny; i++, pc += nx)
          *pc = value;
  }

  void DijkstraPlanner::clearRobotCell(const geometry_msgs::PoseStamped& global_pose, unsigned int mx, unsigned int my)
  {
      if (!initialized_)
      {
          ROS_ERROR("This planner has not been initialized yet, but it is being used, please call initialize() before use");
          return;
      }

      //set the associated costs in the cost map to be free
      costmap_->setCost(mx, my, costmap_2d::FREE_SPACE);
  }

  void DijkstraPlanner::mapToWorld(double mx, double my, double& wx, double& wy)
  {
      wx = costmap_->getOriginX() + (mx+convert_offset_) * costmap_->getResolution();
      wy = costmap_->getOriginY() + (my+convert_offset_) * costmap_->getResolution();
  }

  bool DijkstraPlanner::worldToMap(double wx, double wy, double& mx, double& my)
  {
      double origin_x = costmap_->getOriginX(), origin_y = costmap_->getOriginY();
      double resolution = costmap_->getResolution();

      if (wx < origin_x || wy < origin_y)
          return false;

      mx = (wx - origin_x) / resolution - convert_offset_;
      my = (wy - origin_y) / resolution - convert_offset_;

      if (mx < costmap_->getSizeInCellsX() && my < costmap_->getSizeInCellsY())
          return true;

      return false;
  }

  void DijkstraPlanner::publishPotential(float* potential)
  {
      int nx = costmap_->getSizeInCellsX(), ny = costmap_->getSizeInCellsY();
      double resolution = costmap_->getResolution();
      nav_msgs::OccupancyGrid grid;
      // Publish Whole Grid
      grid.header.frame_id = frame_id_;
      grid.header.stamp = ros::Time::now();
      grid.info.resolution = resolution;

      grid.info.width = nx;
      grid.info.height = ny;

      double wx, wy;
      costmap_->mapToWorld(0, 0, wx, wy);
      grid.info.origin.position.x = wx - resolution / 2;
      grid.info.origin.position.y = wy - resolution / 2;
      grid.info.origin.position.z = 0.0;
      grid.info.origin.orientation.w = 1.0;

      grid.data.resize(nx * ny);

      float max = 0.0;
      for (unsigned int i = 0; i < grid.data.size(); i++)
      {
          float po = potential[i];
          if (po < POT_HIGH)
          {
              if (po > max)
              {
                  max = po;
              }
          }
      }

      for (unsigned int i = 0; i < grid.data.size(); i++)
      {
          if (potential[i] >= POT_HIGH)
          {
              grid.data[i] = -1;
          } else
              grid.data[i] = potential[i] * publish_scale_ / max;
      }
      potential_pub_.publish(grid);
  }

  bool DijkstraPlanner::makePlan(const geometry_msgs::PoseStamped set_start,
                                 const geometry_msgs::PoseStamped set_goal,
                                 std::vector<geometry_msgs::PoseStamped>& plan)
  {
      boost::mutex::scoped_lock lock(mutex_);
      if (!initialized_)
      {
          ROS_ERROR("This planner has not been initialized yet, but it is being used, please call initialize() before use");
          return false;
      }

      geometry_msgs::PoseStamped start_new;
      geometry_msgs::PoseStamped goal_new;

      this->getNearFreePoint(set_start,start_new,this->default_tolerance_);
      this->getNearFreePoint(set_goal,goal_new,this->default_tolerance_);

      //clear the plan, just in case
      plan.clear();

      std::string global_frame = frame_id_;

      //until tf can handle transforming things that are way in the past... we'll require the goal to be in our global frame
      if (goal_new.header.frame_id != global_frame)
      {
          ROS_ERROR("The goal pose passed to this planner must be in the %s frame.  It is instead in the %s frame.", global_frame.c_str(), goal_new.header.frame_id.c_str());
          return false;
      }

      if (start_new.header.frame_id != global_frame)
      {
          ROS_ERROR("The start pose passed to this planner must be in the %s frame.  It is instead in the %s frame.", global_frame.c_str(), start_new.header.frame_id.c_str());
          return false;
      }

      double wx = start_new.pose.position.x;
      double wy = start_new.pose.position.y;

      unsigned int start_x_i, start_y_i, goal_x_i, goal_y_i;
      double start_x, start_y, goal_x, goal_y;

      if (!costmap_->worldToMap(wx, wy, start_x_i, start_y_i))
      {
          ROS_WARN("The robot's start position is off the global costmap. Planning will always fail, are you sure the robot has been properly localized?");
          return false;
      }

      worldToMap(wx, wy, start_x, start_y);

      wx = goal_new.pose.position.x;
      wy = goal_new.pose.position.y;

      if (!costmap_->worldToMap(wx, wy, goal_x_i, goal_y_i))
      {
          ROS_WARN_THROTTLE(1.0,"The goal sent to the global planner is off the global costmap. Planning will always fail to this goal.");
          return false;
      }
      worldToMap(wx, wy, goal_x, goal_y);

      //clear the starting cell within the costmap because we know it can't be an obstacle
      clearRobotCell(start_new, start_x_i, start_y_i);

      int nx = costmap_->getSizeInCellsX(), ny = costmap_->getSizeInCellsY();

      //make sure to resize the underlying array that Navfn uses
      p_calc_->setSize(nx, ny);
      planner_->setSize(nx, ny);
      path_maker_->setSize(nx, ny);
      potential_array_ = new float[nx * ny];

      outlineMap(costmap_->getCharMap(), nx, ny, costmap_2d::LETHAL_OBSTACLE);

      uint8_t plan_timer = 0;
      planner_->setSafetyControl(true);
      while(plan_timer < 2 && plan.empty())
      {
          //计算可行点矩阵，距离机器人越近值越小
          bool found_legal = planner_->calculatePotentials(costmap_->getCharMap(),
                                                           start_x, start_y, goal_x, goal_y,
                                                           nx * ny * 2, potential_array_);

          planner_->clearEndpoint(costmap_->getCharMap(), potential_array_, goal_x_i, goal_y_i, 2);

      //    this->publishPotential(potential_array_);

          if (found_legal)
          {
              //extract the plan
              //从可行点矩阵提取路径
              if (getPlanFromPotential(start_x, start_y, goal_x, goal_y, goal_new, plan))
              {
                  //make sure the goal we push on has the same timestamp as the rest of the plan
                  geometry_msgs::PoseStamped goal_copy = goal_new;
                  goal_copy.header.stamp = ros::Time::now();
                  plan.push_back(goal_copy);
              }
              else
              {
                if(plan_timer == 1)
                  ROS_ERROR("local planner: Failed to get a plan. This shouldn't happen.");
              }
          }
          else
          {
             if(plan_timer == 1)
              ROS_ERROR("local planner:Failed to get a plan.");
          }
          planner_->setSafetyControl(false);
          plan_timer++;
      }
      if(!plan.empty())
      {
        // optimization path
        optimizationPath(plan,M_PI/10);
        // add orientations if needed
        optimizationOrientation(plan);
        //orientation_filter_->processPath(start_new, plan);
        //publish the plan for visualization purposes
        delete potential_array_;
      }
      return !plan.empty();
  }

  int DijkstraPlanner::optimizationPath(std::vector<geometry_msgs::PoseStamped>& plan,double movement_angle_range)
  {
    if(plan.empty())
      return 0;
    size_t pose_size = plan.size() - 1;
    double px,py,cx,cy,nx,ny,a_p,a_n;
    bool is_run = false;
    int ci = 0;
    for(ci=0;ci<1000;ci++)
    {
      is_run = false;
      for(size_t i=1;i<pose_size;i++)
      {
        px = plan[i-1].pose.position.x;
        py = plan[i-1].pose.position.y;

        cx = plan[i].pose.position.x;
        cy = plan[i].pose.position.y;

        nx = plan[i+1].pose.position.x;
        ny = plan[i+1].pose.position.y;

        a_p = normalizeAngle(atan2(cy-py,cx-px),0,2*M_PI);
        a_n = normalizeAngle(atan2(ny-cy,nx-cx),0,2*M_PI);

        if(std::max(a_p,a_n)-std::min(a_p,a_n) > movement_angle_range)
        {
          plan[i].pose.position.x = (px + nx)/2;
          plan[i].pose.position.y = (py + ny)/2;
          is_run = true;
        }
      }
      if(!is_run)
        return ci;
    }
    return ci;
  }

  void DijkstraPlanner::optimizationOrientation(std::vector<geometry_msgs::PoseStamped> &plan)
  {
    size_t num = plan.size()-1;
    if(num < 1)
      return;
    for(size_t i=0;i<num;i++)
    {
      plan[i].pose.orientation = tf::createQuaternionMsgFromYaw( atan2( plan[i+1].pose.position.y - plan[i].pose.position.y,
                                                                 plan[i+1].pose.position.x - plan[i].pose.position.x ) );
    }
  }

  double DijkstraPlanner::distance(double x1,double y1,double x2,double y2)
  {
    return sqrt( (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2) );
  }

  bool DijkstraPlanner::isAroundFree(unsigned int mx, unsigned int my)
  {
    if(mx <= 1 || my <= 1 || mx >= this->costmap_->getSizeInCellsX()-1 || my >= this->costmap_->getSizeInCellsY()-1)
      return false;
    int x,y;
    for(int i=-1;i<=1;i++)
    {
      for(int j=-1;j<=1;j++)
      {
        x = static_cast<int>(mx) + i;
        y = static_cast<int>(my) + j;
        if(this->costmap_->getCost(static_cast<unsigned int>(x),static_cast<unsigned int>(y)) != costmap_2d::FREE_SPACE)
          return false;
      }
    }
    return true;
  }

  void DijkstraPlanner::getNearFreePoint(const geometry_msgs::PoseStamped in,
                                       geometry_msgs::PoseStamped& out,
                                       double tolerance)
  {
    out = in;
    unsigned int grid_size = static_cast<unsigned int>(tolerance/costmap_->getResolution() + 0.5);
    if(grid_size<1)
    {
      out = in;
      return;
    }

    unsigned int mx0,my0;
    if(costmap_->worldToMap(in.pose.position.x,in.pose.position.y,mx0,my0))
    {
      if(this->isAroundFree(mx0,my0))
        return;
      unsigned int minx,maxx,miny,maxy;
      double wx = 0.0,wy = 0.0;
      double min_move_cost = 10000000.0;
      minx = mx0-grid_size>0?mx0-grid_size:0;
      maxx = mx0+grid_size<costmap_->getSizeInCellsX()?mx0+grid_size:costmap_->getSizeInCellsX();
      miny = my0-grid_size>0?my0-grid_size:0;
      maxy = my0+grid_size<costmap_->getSizeInCellsY()?my0+grid_size:costmap_->getSizeInCellsY();
      for(unsigned int i=minx;i<=maxx;i++)
      {
        for(unsigned int j=miny;j<=maxy;j++)
        {
          costmap_->mapToWorld(i,j,wx,wy);
          double current_move_cost = this->distance(in.pose.position.x,in.pose.position.y,wx,wy);
          if(!this->isAroundFree(i,j) || current_move_cost > tolerance)
            continue;
          if(min_move_cost > current_move_cost)
          {
            min_move_cost = current_move_cost;
            out.pose.position.x = wx;
            out.pose.position.y = wy;
          }
        }
      }
    }
  }

  bool DijkstraPlanner::getPlanFromPotential(double start_x, double start_y,
                                           double goal_x, double goal_y,
                                           const geometry_msgs::PoseStamped& goal,
                                           std::vector<geometry_msgs::PoseStamped>& plan)
  {
      if (!initialized_)
      {
          ROS_ERROR("This planner has not been initialized yet, but it is being used, please call initialize() before use");
          return false;
      }

      //clear the plan, just in case
      plan.clear();

      std::vector<std::pair<float, float> > path;

      if (!path_maker_->getPath(potential_array_, start_x, start_y, goal_x, goal_y, path))
      {
          ROS_ERROR("NO PATH!");
          return false;
      }

      ros::Time plan_time = ros::Time::now();
      int path_size_num = path.size() -1;
      for (int i = path_size_num; i>=0; i--)
      {
          std::pair<float, float> point = path[i];
          //convert the plan to world coordinates
          double world_x, world_y;
          mapToWorld(point.first, point.second, world_x, world_y);

          geometry_msgs::PoseStamped pose;
          pose.header.stamp = plan_time;
          pose.header.frame_id = frame_id_;
          pose.pose.position.x = world_x;
          pose.pose.position.y = world_y;
          pose.pose.position.z = 0.0;
          pose.pose.orientation.x = 0.0;
          pose.pose.orientation.y = 0.0;
          pose.pose.orientation.z = 0.0;
          pose.pose.orientation.w = 1.0;

          if(i != path_size_num)
          {
            double cx,cy,px,py;
            cx = pose.pose.position.x;
            cy = pose.pose.position.y;
            px = plan.back().pose.position.x;
            py = plan.back().pose.position.y;
            if( sqrt( (cx-px)*(cx-px) + (cy-py)*(cy-py) ) > 0.05)
            {
              geometry_msgs::PoseStamped pose_insert = pose;
              pose_insert.pose.position.x = (cx+px)/2;
              pose_insert.pose.position.y = (cy+py)/2;
              plan.push_back(pose_insert);
            }
          }
          plan.push_back(pose);
      }
      return !plan.empty();
  }


};
