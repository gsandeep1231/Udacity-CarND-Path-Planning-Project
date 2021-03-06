# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

### Project Details
#### 1. Get vehicle reference location
Depending on whether the previous path of the vehicle has any points or not, we initialize the reference location of the vehicle. If it does, we choose the last two point from previous path and use them, else use the vehicle's location from the senor fusion data and extrapolate its location in previous time period.
We also initialize two values for our vectors ptsx and ptsy that will be later used to create spline.

			double ref_x = car_x;
			double ref_y = car_y;
			double ref_yaw = car_yaw;
			
			if (path_size < 2) {
				ptsx.push_back(car_x - cos(car_yaw));
				ptsx.push_back(car_x);
				
				ptsy.push_back(car_y - sin(car_yaw));
				ptsy.push_back(car_y);
			} else {
				ref_x = previous_path_x[path_size-1];
				ref_y = previous_path_y[path_size-1];
				double prev_ref_x = previous_path_x[path_size-2];
				double prev_ref_y = previous_path_y[path_size-2];
				ref_yaw = atan2(ref_y-prev_ref_y, ref_x-prev_ref_x);
				
				ptsx.push_back(prev_ref_x);
				ptsx.push_back(ref_x);
				
				ptsy.push_back(prev_ref_y);
				ptsy.push_back(ref_y);
			}

#### 2. Get more waypoints
Using frenet coordinates, estimate vehicles's future waypoints by only incrementing s and using map waypoints.
Convert this to global coordinates using inbuilt function getXY.
After trying different values of incremental s for future waypoints, I found that adding 50 to the vehicle reference position gave better jerk acceleration.

			vector<double> next_wp0 = getXY(car_s+50,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp1 = getXY(car_s+100,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp2 = getXY(car_s+120,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

#### 3. Shift car reference
Shifting the car reference to zero makes the calculations easier to estimate trajectory.

			for(int i=0; i<ptsx.size(); i++) {
				double shift_x = ptsx[i] - ref_x;
				double shift_y = ptsy[i] - ref_y;
				
				ptsx[i] = (shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw));
				ptsy[i] = (shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw));
			}
			
#### 4. Create a spline
Using spline.h from the project description, create a spline using the vectors ptsx and ptsy.

			tk::spline s;
			s.set_points(ptsx, ptsy);
			
#### 5. Create next trajectory values
Finally using all the points from previous trajectory, initialize the current trajectory.

			vector<double> next_x_vals;
			vector<double> next_y_vals;
          	
			for(int i=0; i<path_size; i++) {
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}

Since we expect 50 points in the trajectory, for the remaining ones estimate the trajectory points based on the reference velocity.
Using the reference velocity find the number of intervals (N) needed to cover the distance if each interval is 20ms apart. Finally using this N calculate the next x and y location of the vehicle. Since we earlier moved the reference to zero, revert the reference back to original location to get the actual x and y location of the vehicle in terms of global coordinates.

			double target_x = 30.0;
			double target_y = s(target_x);
			double target_dist = sqrt((target_x*target_x) + (target_y*target_y));
			
			double x_add_on = 0;
			for (int i=0; i<50-path_size; i++) {
				double N = target_dist/(0.02*ref_vel/2.24);
				double x_point = x_add_on+(target_x/N);
				double y_point = s(x_point);
				x_add_on = x_point;
				
				double x_ref = x_point;
				double y_ref = y_point;
				
				// rotate back to normal
				x_point = (x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw));
				y_point = (x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw));
				
				x_point += ref_x;
				y_point += ref_y;
				
				next_x_vals.push_back(x_point);
				next_y_vals.push_back(y_point);
			}
#### 6. Selecting reference velocity and lane
Reference velocity is initially set to zero. When the vehicle start, the velocity is slowly incremented until it reaches the speed limit of 49.5mph.
We use three flags to determine how our model behaves and changes the reference velocity and lane of our vehicle using data from sensor fusion.

too_close:
If any vehicle in front of our vehicle, in the same lane, has its end path within 30m of our end path, then we set this flag to 1. This then indicates that the vehicle in front of us is too close.

				if (d < (2+4*lane+2) && d > (2+4*lane-2)) {
					if((other_car_end_s > car_end_s) && ((other_car_end_s-car_end_s) < 30)) {
						front_car_speed = other_car_speed;
						too_close = true;
					}
				}

safe_left:
This is intialized to 1. If we find that we are in the leftmost lane or if any vehicle is found to our left then we set this flag to 0.

				if (lane>0) {
					int new_lane = lane-1;
					if (d < (2+4*new_lane+2) && d > (2+4*new_lane-2)) {
						cout << "vehicle found with id: " << sensor_fusion[i][0] << endl;
						if ((other_car_s < car_s && (other_car_end_s) > car_s) || (other_car_s > car_s && other_car_s < (car_end_s+10))) {
							safe_left = false;
						}
					}
				} else {
					safe_left = false;
				}
				
safe_right:
This is intialized to 1. If we find that we are in the rightmost lane or if any vehicle is found to our right then we set this flag to 0.

				if (lane<2) {
					int new_lane = lane+1;
					if (d < (2+4*new_lane+2) && d > (2+4*new_lane-2)) {
						cout << "vehicle found with id: " << sensor_fusion[i][0] << endl;
						if ((other_car_s < car_s && (other_car_end_s) > car_s) || (other_car_s > car_s && other_car_s < (car_end_s+10))) {
							safe_right = false;
						}
					}
				} else {
					safe_right = false;
				}

Using these flags, I then create a FSM model that has three states, keep lane, move left or move right.
This is depicted in below conditional statements by modifying the reference velocity and lane number.

			if (too_close) {
				ref_vel -= 0.224;
				if (ref_vel < front_car_speed) {
					ref_vel = front_car_speed;
				}
				if (safe_left) {
					lane--;
					wait_time = 0;
				} else if (safe_right && wait_time > 100) {
					lane++;
					wait_time = 0;
				} else {
					wait_time++;
				}
			} else if (ref_vel<49.5){
				ref_vel += 0.224;
			}

### Simulator.
You can download the Term3 Simulator which contains the Path Planning Project from the [releases tab (https://github.com/udacity/self-driving-car-sim/releases).

### Goals
In this project your goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. You will be provided the car's localization and sensor fusion data, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 50 m/s^3.

#### The map of the highway is in data/highway_map.txt
Each waypoint in the list contains  [x,y,s,dx,dy] values. x and y are the waypoint's map coordinate position, the s value is the distance along the road to get to that waypoint in meters, the dx and dy values define the unit normal vector pointing outward of the highway loop.

The highway's waypoints loop around so the frenet s value, distance along the road, goes from 0 to 6945.554.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.

Here is the data provided from the Simulator to the C++ Program

#### Main car's localization Data (No Noise)

["x"] The car's x position in map coordinates

["y"] The car's y position in map coordinates

["s"] The car's s position in frenet coordinates

["d"] The car's d position in frenet coordinates

["yaw"] The car's yaw angle in the map

["speed"] The car's speed in MPH

#### Previous path data given to the Planner

//Note: Return the previous list but with processed points removed, can be a nice tool to show how far along
the path has processed since last time. 

["previous_path_x"] The previous list of x points previously given to the simulator

["previous_path_y"] The previous list of y points previously given to the simulator

#### Previous path's end s and d values 

["end_path_s"] The previous list's last point's frenet s value

["end_path_d"] The previous list's last point's frenet d value

#### Sensor Fusion Data, a list of all other car's attributes on the same side of the road. (No Noise)

["sensor_fusion"] A 2d vector of cars and then that car's [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates. 

## Details

1. The car uses a perfect controller and will visit every (x,y) point it recieves in the list every .02 seconds. The units for the (x,y) points are in meters and the spacing of the points determines the speed of the car. The vector going from a point to the next point in the list dictates the angle of the car. Acceleration both in the tangential and normal directions is measured along with the jerk, the rate of change of total Acceleration. The (x,y) point paths that the planner recieves should not have a total acceleration that goes over 10 m/s^2, also the jerk should not go over 50 m/s^3. (NOTE: As this is BETA, these requirements might change. Also currently jerk is over a .02 second interval, it would probably be better to average total acceleration over 1 second and measure jerk from that.

2. There will be some latency between the simulator running and the path planner returning a path, with optimized code usually its not very long maybe just 1-3 time steps. During this delay the simulator will continue using points that it was last given, because of this its a good idea to store the last points you have used so you can have a smooth transition. previous_path_x, and previous_path_y can be helpful for this transition since they show the last points given to the simulator controller with the processed points already removed. You would either return a path that extends this previous path or make sure to create a new path that has a smooth transition with this last path.

## Tips

A really helpful resource for doing this project and creating smooth trajectories was using http://kluge.in-chemnitz.de/opensource/spline/, the spline function is in a single hearder file is really easy to use.

---

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).

## Project Instructions and Rubric

Note: regardless of the changes you make, your project must be buildable using
cmake and make!


## Call for IDE Profiles Pull Requests

Help your fellow students!

We decided to create Makefiles with cmake to keep this project as platform
agnostic as possible. Similarly, we omitted IDE profiles in order to ensure
that students don't feel pressured to use one IDE or another.

However! I'd love to help people get up and running with their IDEs of choice.
If you've created a profile for an IDE that you think other students would
appreciate, we'd love to have you add the requisite profile files and
instructions to ide_profiles/. For example if you wanted to add a VS Code
profile, you'd add:

* /ide_profiles/vscode/.vscode
* /ide_profiles/vscode/README.md

The README should explain what the profile does, how to take advantage of it,
and how to install it.

Frankly, I've never been involved in a project with multiple IDE profiles
before. I believe the best way to handle this would be to keep them out of the
repo root to avoid clutter. My expectation is that most profiles will include
instructions to copy files to a new location to get picked up by the IDE, but
that's just a guess.

One last note here: regardless of the IDE used, every submitted project must
still be compilable with cmake and make./

## How to write a README
A well written README file can enhance your project and portfolio.  Develop your abilities to create professional README files by completing [this free course](https://www.udacity.com/course/writing-readmes--ud777).

