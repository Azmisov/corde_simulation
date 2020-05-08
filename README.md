# CoRdE Implementation
This is an implementation of the paper "CoRdE: Cosserat Rod Elements for the Dynamic Simulation of One-Dimensional Elastic Objects" (2007, J. Spillman, M. Teschner). It can be used to simulate things like rope, string, chains, elastic rods, wires, bars, springs, and more. I originally intended to use it for 3D rope simulation in games, but could be adapted to these other use cases. I am only including the header file `Rope.hpp`, which contains all the relevant simulation code. It is using [xtensor](https://github.com/xtensor-stack/xtensor) for matrix operations (with my own aliases for `xm::vec3`, `xm::mat4` etc contained in `xutils.hpp`).

![Screencast](https://raw.githubusercontent.com/Azmisov/corde_simulation/master/screencast.gif)

Here's how it could be used, once included in your own project:
```c++
// constructor automatically generates a simple ine
Rope r;
// add external forces/torques
// note hwoever that currently I have gravity + ground collision + damping + air drag hard coded in the code
r.pts.back().F += (xm::vec3) {1,0,0};
// this performs a simulation timestep
r.update();
// there are a couple builtin vertex array objects for vizualizing
// this updates the vao's with new vertex/line data
r.update_buffers();
// here's what drawing might look like:
uint32_t rpts = rope.pts.size();
vertex_only_shader.use();
glUniform3f(uniform_color_attribute,1.f,.5f,0.f);
glBindVertexArray(rope.segments_vao);
glDrawArrays(GL_LINE_STRIP, 0, rpts);
// rope points
glUniform3f(uniform_color_attribute,1.f,.85f,0.f);
glDrawArrays(GL_POINTS, 0, rpts);
// rope directors
vertex_and_color_shader.use();
glBindVertexArray(rope.directors_vao);
glDrawArrays(GL_LINES, 0, (rpts-1)*6);
```

## Implementation Notes
I'm mostly uploading this in case it is helpful to someone else. I could not find any existing public implementations myself. However, this simulation method is not without problems. First of all, the paper leaves out quite a few details, so I had to guess on some of the mathematics. The paper was written a while ago, so the authors couldn't provide much assistance either. I eventually was able to get the simulation working and seems to be somewhat stable. See the `notes` folder for:

- A copy of the CoRdE technical paper
- The author's dissertation, which includes a few additional clarifying details
- My own notes on the paper. I documented most of my guesses about the mathematics here, as well as a number of problems I encountered with the method. For instance, the material constants did not give good results as given in the paper; I had to play with them a little, and they still are not very intuitive. Besides a few smaller mathematical details that I was not sure of (especially handling the quaternions), the biggest change I made was in the partial derivative calculations. One must take the partial derivatives of the some energy equations to compute internal forces and torques for the rope. However, the author's do this using automatic differentiation in maple, and don't provide thorough details. I used Sage to compute these derivatives, and then reduce the resulting equations to simplified forms myself manually. I believe my manual derivations for the partial derivatives are much faster to compute.
	- `rope_math.ipynb`: the Jupyter Sage notebook I used to aid in the partial derivatives. I also verify that my manual simplifications are still correct here.
	- `[name]_pd.txt`: I used these files to record my simplifications of the Sage output expressions

**Now, I don't really advocate using this method for actual simulation.** This repo is more of a reference for anyone having trouble with similar simulations or just for the curious soul. Here is why:
- The simulation method is _explicit_, meaning it can need very small timesteps for the simulation to be stable, especially if you are simulating stiffer objects. I believe my implementation is quite efficient though, so its pretty good with say 50 segments. If you put it on the GPU (which is surely possible with this method), you could do much more. However, there are newer methods that are faster now.
- The material constants aren't intuitive. I couldn't get them to match up with what is listed in the original paper though, so if you found a way to fix this, perhaps this wouldn't be a problem.
- Perhaps due to timestep or improper material constants, I found it still became unstable in many situations. More fine tuning and work would be needed to figure out where this is coming from and possibly fix any errors in the code.
- Its an old method (2007), and the original authors likely cannot provide help if you are having problems
- Collision handling and constraint handling are apparently not quite as easy as other competing methods (e.g. position based)
- There are better methods available, at least for realtime simulation (as was my use case), that have better integration with other simulation use caes. I recommend checking out these methods instead: "Position and Orientation Based Cosserat Rods" (2016), "Realtime simulation of inextensible surgical thread using a Kirchhoff rod model" (2017), and "Direct Position-Based Solver for Stiff Rods" (2018). There may be newer papers available as well. I also recently discovered the open source [PositionBasedDynamics](https://github.com/InteractiveComputerGraphics/PositionBasedDynamics) project, which implements several of these papers; go check it out.
