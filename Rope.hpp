#pragma once
#include <vector>
#include <cmath>
#include "xutils.hpp"

class Rope;

/* Will be used to store mass points + directors together. The last point's director (and related vars)
	will be a dummy (since directors are in between two mass points), but this way I think
	will be cleaner than having the two separated. Also will be cleaner for extending the rope
*/
struct RopePoint{
	xm::vec3 r, // position
		v, // translational velocity
		w; // rotational velocity
	xm::vec4 q; // director quaternion
	// material constants
	float m, // mass (kg)
		li, // resting length
		lj; // resting torsion

	// given F/T etc accumulated from simulation, updates current state
	void update(const Rope &rope, float dt, bool islast, uint32_t idx);

	// These are only needed on a rolling window as we compute internal forces
	xm::vec3 F; // forces
	xm::vec4 T; // torques
	xm::vec3 Bwq2; // for diss torque calculation
	void update_temps(){
		// Bwq2 = B0q*Q*w
		// would *.5, but we're absorbing it into the larger calculation
		Bwq2 = xt::linalg::dot(B0q(),xt::linalg::dot(Qmat(),w));
	}

	// quaternion matrix helpers; a few of these are magic matrices found from the partial derivatives
	// matrix of B0*q
	inline xm::mat34 B0q(){
		return {
			{ q[3],-q[2], q[1],-q[0]}, // B0_1
			{ q[2], q[3],-q[0],-q[1]}, // B0_2
			{-q[1], q[0], q[3],-q[2]} // B0_3
		};
	}
	// quaternion matrix w/ first column chopped off
	// TODO: verify that one entry would have been zero
	inline xm::mat43 Qmat(){
		return {
			// standard wiki matrix
			//{-q[1],-q[2],-q[3]},
			//{ q[0],-q[3], q[2]},
			//{ q[3], q[0],-q[1]},
			//{-q[2], q[1], q[0]}

			// w/ first column chopped off
			//{ q[2],-q[1], q[0]},
			//{ q[3], q[0], q[1]},
			//{-q[0], q[3], q[2]},
			//{-q[1],-q[2], q[3]}
			
			// this is params given in dissertaiton
			//{ q[3], q[2],-q[1], q[0]},
			//{-q[2], q[3], q[0], q[1]},
			//{ q[1],-q[0], q[3], q[2]},
			//{-q[0],-q[1],-q[2], q[3]}

			// based on dissertation, chopping off last column
			{ q[3], q[2],-q[1]},
			{-q[2], q[3], q[0]},
			{ q[1],-q[0], q[3]},
			{-q[0],-q[1],-q[2]}
		};
	}
	inline xm::mat4 Qmat_full(){
		return {
			{ q[3], q[2],-q[1], q[0]},
			{-q[2], q[3], q[0], q[1]},
			{ q[1],-q[0], q[3], q[2]},
			{-q[0],-q[1],-q[2], q[3]}
		};
	}
	// swizzle for bending torque
	inline xm::mat43 Qmat_bend(){
		return {
			{-q[3], q[2],-q[1]},
			{-q[2],-q[3], q[0]},
			{ q[1],-q[0],-q[3]},
			{ q[0], q[1], q[2]}
		};
	}
	// swizzle for constraint torque
	inline xm::mat43 Qmat_cons(){
		return {
			{-q[2], q[3], q[0]},
			{-q[3],-q[2], q[1]},
			{-q[0],-q[1],-q[2]},
			{-q[1], q[0],-q[3]}
		};
	}

	// converts quaternion to directors (just rotation matrix)
	xm::mat3 directors(){
		float aa = q[0]*q[0], ab = q[0]*q[1], ac = q[0]*q[2], ad = q[0]*q[3],
			bb = q[1]*q[1], bc = q[1]*q[2], bd = q[1]*q[3],
			cc = q[2]*q[2], cd = q[2]*q[3];
		return {	
			{1-2*(bb+cc), 2*(ab-cd), 2*(ac+bd)},
			{2*(ab+cd), 1-2*(aa+cc), 2*(bc-ad)},
			{2*(ac-bd), 2*(bc+ad), 1-2*(aa+bb)}
		};
	}
	// gets quaternion, given a euclidean R3 orthonormal basis
	static xm::vec4 eucToQ(const xm::mat3 &em){
		//*
		xm::vec4 qp;
		float t;
		if (em(2,2) < 0){
			if (em(0,0) > em(1,1)){
				t = 1+em(0,0)-em(1,1)-em(2,2);
				qp = {t, em(0,1)+em(1,0), em(2,0)+em(0,2), em(2,1)-em(1,2)};
			}
			else{
				t = 1-em(0,0)+em(1,1)-em(2,2);
				qp = {em(0,1)+em(1,0), t, em(1,2)+em(2,1), em(0,2)-em(2,0)};
			}
		}
		else{
			if (em(0,0) < -em(1,1)){
				t = 1-em(0,0)-em(1,1)+em(2,2);
				qp = {em(2,0)+em(0,2), em(1,2)+em(2,1), t, em(1,0)-em(0,1)};
			}
			else{
				t = 1+em(0,0)+em(1,1)+em(2,2);
				qp = {em(1,2)-em(2,1), em(2,0)-em(0,2), em(0,1)-em(1,0), -t};
			}
		}
		qp *= 0.5/std::sqrt(t);
		return qp;
		//*/
	}
};

class Rope{
public:
	uint32_t segments_vao, segments_vbo,
		directors_vao, directors_vbo;
	// material properties
	float r = .01, // radius (m)
		rho = 1300, // density (kg/m^3)
		E = .5*1e3, // young modulus (MPa), new gpa
		G = .5*1e3, // shearing modulus (MPa)
		Es = 20*1e6, // stretch modulus (MPa)
		k = 100*1e3, // spring constant (10^3*kg*m/s^2)
		gt = 10*1e-6, // translation internal friction constant (10^-6*kg*m^3/s)
		gr = 1*1e-6; //rotation internal friction constant (10^-6*kg*m^3/s)
	// derived properties
	float Ix, Iy, // inertia tensor
		Kx, Ky, // stiffness tensor
		Ks; // stretching stiffness

	std::vector<RopePoint> pts;

	Rope(){
		init_material();
		generate_line(3, 50, {0,-1.5,1.2}, {0,1,0});
	}

	void init_material(){
		float half_pir2 = r*r*1.570796327;
		Iy = rho*half_pir2;
		Ix = Iy/2.;
		Ky = G*half_pir2;
		Kx = Ky/2.;
		Ks = Es*2*half_pir2;
	}

	/* generates a straight line of rope points and initializes them to static initial state
		length: length of line in meters
		divs: how many segments there should be (+1 for points)
		start: starting point of line
		dir: start+dir*length = endpoint
	*/
	void generate_line(float length, uint32_t divs, xm::vec3 start, xm::vec3 dir){
		float v = r*r*3.141592654*length, // volume
			mi = v*rho/divs, // segment mass
			li = length/divs; // segment len

		dir /= xt::norm_l2(dir)();
		xm::vec3 end = start+length*dir;
		
		/* Compute directors
			1st basis = line direction
			we can have 2nd basis vector be fixed to xy plane (z=0), forcing it all to be flat and untwisted:
				xd*x2 + yd*y2 = 0
				y2 = -x2*xd/yd
			now normalize to unit length:
				n = sqrt(x2^2 + x2^2(xd/yd)^2)
					= x2*sqrt(1 + xd*xd/(yd*yd))
				x2 = 1/sqrt(1 + (xd/yd)^2)
				y2 = -xd/(yd*sqrt(1 + (xd/y2)^2))
					= -xd/(sqrt(yd^2 + xd^2))
				the most numerically stable? probably not
			3rd basis can be found from cross product; we'll do right-handed, so if z is negative, we negate 2nd and 3rd
		*/
		xm::mat3 directors;
		auto d1 = xt::view(directors, xt::all(), 0);
		auto d2 = xt::view(directors, xt::all(), 1);
		auto d3 = xt::view(directors, xt::all(), 2);
		d3 = dir; // d3 should match curvature
		if (dir[1] < 1e-12)
			d1 = (xm::vec3) {0,1,0};
		else{
			float div = -dir[0]/dir[1],
				norm = 1/std::sqrt(1+div*div);
			d1 = (xm::vec3) {norm,div*norm,0};
		}
		d2 = xt::linalg::cross(d3,d1);
		if (d2[2] < 0){
			d1 = -d1;
			d2 = -d2;
		}
		xm::vec4 qi = RopePoint::eucToQ(directors);
		
		for (uint32_t i=0; i<=divs; i++){
			pts.emplace_back();
			RopePoint &p = pts.back();
			float t = i/static_cast<float>(divs);
			p.r = start*(1-t) + end*t;
			// start/end get half a segment
			p.m = !i || i==divs ? .5*mi : mi;
			// at rest initially
			p.v.fill(0);
			p.w.fill(0);
			p.F.fill(0);
			p.T.fill(0);
			// these are the centerline props
			if (i != divs){
				p.q = qi;
				// for arbitrary points, can calculate as:
				//p.li = xt::eval(xt::norm_l2(p.r - other.r)())
				// since all even segments, and no initial bending, lj=li
				p.li = p.lj = li;
			}
		}

		// OpenGL initailization
		// segments
		glGenVertexArrays(1, &segments_vao);
		glBindVertexArray(segments_vao);
		glGenBuffers(1, &segments_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, segments_vbo);
		glBufferData(GL_ARRAY_BUFFER, pts.size()*3*sizeof(float), nullptr, GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0*sizeof(float), (void*)(0*sizeof(float))); // pos
		glEnableVertexAttribArray(0);

		// directors
		uint32_t dpts = (pts.size()-1)*36; // 3 lines * 2 points * (3 axes + 3 cols)
		glGenVertexArrays(1, &directors_vao);
		glBindVertexArray(directors_vao);
		glGenBuffers(1, &directors_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, directors_vbo);
		glBufferData(GL_ARRAY_BUFFER, dpts*sizeof(float), nullptr, GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)(0*sizeof(float))); // pos
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)(3*sizeof(float))); // col
		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);
	}
	// fills vbo buffers with updated data
	void update_buffers(){
		// segments (and points)
		glBindBuffer(GL_ARRAY_BUFFER, segments_vbo);
		float *sb = static_cast<float*>(glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY));
		uint32_t si = 0;
		for (auto &p: pts){
			for (uint8_t j=0; j<3; j++, si++){
				sb[si] = p.r[j];
			}
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);

		// directors; rgb axis vizualization
		glBindBuffer(GL_ARRAY_BUFFER, directors_vbo);
		float *db = static_cast<float*>(glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY));
		// how long the director axis viz should be, relative to segment
		const float dir_prop = .25;
		uint32_t di = 0;
		for (uint32_t pi=0; pi<pts.size()-1; pi++){
			auto &a = pts[pi], &b = pts[pi+1];
			float len = a.li*dir_prop;
			xm::vec3 pos = .5*(a.r+b.r);
			xm::mat3 dirs = a.directors();
			// create viz for each axis
			for (uint8_t ax=0; ax<3; ax++){
				for (uint8_t j=0; j<3; j++, di++){
					float c = j == ax ? 1 : 0;
					// start point
					db[di] = pos[j];
					db[di+3] = c;
					// end point
					db[di+6] = dirs(j,ax)*len + pos[j];
					db[di+9] = c;
				}
				di += 9;
			}
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);
	}

	// This is where simulation happens
	void update(){
		float dt = .1/1000.; // .1ms

		// TODO, maybe some assertions that sqrts are  > 0
		pts[0].update_temps();
		uint32_t last_segment = pts.size()-2;
		for (size_t i=0; i<=last_segment; i++){
			RopePoint &a = pts[i], &b = pts[i+1];
			b.update_temps();

			// Components
			const xm::vec3 dr = b.r - a.r, dv = b.v - a.v;
			const xm::vec4 dq = b.q - a.q, sq = b.q + a.q;
			const float drr = xt::norm_sq(dr)(),
				dr_l2_inv = 1/std::sqrt(drr),
				drv = xt::linalg::dot(dr,dv)(),
				li5_inv = std::pow(a.li,-5);

			// Compute internal forces/torques
			// TODO: Is there a relationionship between B, Q2, d3, Q? Maybe Q1 as well?
			
			// stretch force
			xm::vec3 F_stretch = Ks*(dr_l2_inv - 1/a.li)*dr; // TODO: precompute ks*(1-1/li) ?
			assert(!xt::any(xt::isnan(F_stretch)) && "stretch nan");
			a.F -= F_stretch;
			b.F += F_stretch;
			
			// bending torque
			if (i != last_segment){
				// matrix of B*(qi+qj)
				const xm::mat34 B = {
					{ sq[3], sq[2],-sq[1],-sq[0]},
					{-sq[2], sq[3], sq[0],-sq[1]},
					{ sq[1],-sq[0], sq[3],-sq[2]}
				};
				// assuming resting bend, uk = 0
				xm::vec3 K = {Kx,Kx,Ky};
				xm::vec3 bend_comp = 2/a.lj*K*xt::linalg::dot(B,dq);
				assert(!xt::any(xt::isnan(bend_comp)) && "bend nan");
				a.T -= xt::linalg::dot(b.Qmat_bend(), bend_comp);
				b.T += xt::linalg::dot(a.Qmat_bend(), bend_comp);
			}
			
			// d3 constraint force
			xm::vec3 d3 = {
				2*(a.q[0]*a.q[2] + a.q[1]*a.q[3]),
				2*(a.q[1]*a.q[2] - a.q[0]*a.q[3]),
				//a.q[2]*a.q[2] + a.q[3]*a.q[3] - a.q[0]*a.q[0] - a.q[1]*a.q[1]
				// this is the calculation without substituting norm(q) for 1:
				1-2*(a.q[0]*a.q[0] + a.q[1]*a.q[1])
			};
			xm::vec3 d3_err = dr*dr_l2_inv - d3;
			xm::vec3 F_constraint = k*a.li*dr_l2_inv*(dr*xt::linalg::dot(dr,d3_err)/drr - d3_err);
			assert(!xt::any(xt::isnan(F_constraint)) && "d3 f nan");
			a.F -= F_constraint;
			b.F += F_constraint;

			// d3 constraint torque
			auto d3_t = 2*k*a.li*xt::linalg::dot(a.Qmat_cons(),d3_err);
			assert(!xt::any(xt::isnan(d3_t)) && "d3 t nan");
			a.T -= d3_t;

			// dissipation force
			xm::vec3 F_diss = -gt*li5_inv*drr*drv*dr; // TODO: maybe combine li5_inv w gt
			assert(!xt::any(xt::isnan(F_diss)) && "f diss nan");
			a.F -= F_diss;
			b.F += F_diss;

			// dissipation torque
			if (i != last_segment){
				// TODO: an idea is to combine the i and i+1 parts into one, then do a single B0t_a multiplication on that
				// TODO: another idea is just splitting out these pde formulas so that we can reuse more of the computations across iterations
				xm::vec3 tdiss_comp = 2*gr/a.lj*(a.Bwq2 - b.Bwq2);
				assert(!xt::any(xt::isnan(tdiss_comp)) && "f diss nan");
				a.T -= xt::linalg::dot(xt::transpose(a.B0q()), tdiss_comp);
				b.T += xt::linalg::dot(xt::transpose(b.B0q()), tdiss_comp);
			}

			// Air drag
			// maybe my calculations are wrong, but seems to oscillate a lot even with air drag; hence extra_drag multiplier
			float extra_drag = 12;
			const float drag_const = -.5*1.2*1.2*r*extra_drag;
			xm::vec3 va = .75*a.v+.25*b.v,
				vb = .25*a.v+.75*b.v;
			float len_va = xt::norm_sq(va)(),
				len_vb = xt::norm_sq(vb)();
			if (len_va > 1e-5){
				va /= std::sqrt(len_va);
				float proj = xt::linalg::dot(dr, va)();				
				proj *= proj;
				if (proj <= drr){
					xm::vec3 drag = drag_const*std::sqrt(drr - proj)*len_va*va;
					assert(!xt::any(xt::isnan(drag)) && "a drag nan");
					a.F += drag;
				}
			}
			if (len_vb > 1e-5){
				vb /= std::sqrt(len_vb);
				float proj = xt::linalg::dot(dr, vb)();
				proj *= proj;
				if (proj <= drr){
					xm::vec3 drag = drag_const*std::sqrt(drr - proj)*len_vb*vb;
					assert(!xt::any(xt::isnan(drag)) && "b drag nan");
					b.F += drag;
				}
			}
		}
		
		// update after; currently using things like p/v for next point, so can't update until after
		for (uint32_t i=0; i<pts.size(); i++)
			pts[i].update(*this, dt, i > last_segment, i);
	}
};

void RopePoint::update(const Rope &rope, float dt, bool islast, uint32_t idx){
	// immovable constraint an first point
	if (idx){
		xm::vec3 gravity = (xm::vec3) {0,0,-9.8*m};
		F += gravity;
		// collision w/ ground
		if (r[2] <= 0){
			r[2] = 0;
			F[2] = 0;
			// kinetic friction
			float kf = .3*gravity[2];
			float slide_norm = std::sqrt(v[0]*v[0]+v[1]*v[1]);
			if (slide_norm > 1e-5){
				F[0] += kf*v[0]/slide_norm;
				F[1] += kf*v[1]/slide_norm;
			}
			// bounce restitution
			v[2] = -v[2]*.5;
			// damping, to sabilize it a bit
			float damp = .93;
			v[0] *= damp;
			v[1] *= damp;
			if (v[2] < 1e-4)
				v[2] = 0;
		}
		v += dt/m*F; // TODO: store inv mass?
		r += dt*v;
		// unilateral damping
	}

	if (!islast){
		// euclidean torques
		// 2 or .5?
		// we want the last entry to be zero, but for some reason it never is
		xm::vec4 t_eucl = 2*xt::linalg::dot(xt::transpose(Qmat_full()), T);

		// torque transducer
		if (!idx){
			float v = 1.;
			t_eucl[2] += v;
		}
		
		// I is diagonal, so equiv to array mult/divide
		//w += dt*(t_eucl - xt::linalg::cross(w, w*rope.I))/rope.I;
		// manual unroll
		// TODO: could store 1/Ix, Iy-Ix, 1/Iy
		// dissertation says multiply by lj
		float aIx = rope.Ix*lj;
		//float dI = w[2]*(rope.Iy - aIx);
		//w[0] += dt*(t_eucl[0] - w[1]*dI)/aIx;
		//w[1] += dt*(t_eucl[1] + w[0]*dI)/aIx;
		//w[2] += dt*t_eucl[2]/rope.Iy;
		// reduction given in dissertaiton
		w[0] += dt*t_eucl[0]/aIx;
		w[1] += dt*t_eucl[1]/aIx;
		w[2] += dt*t_eucl[2]/rope.Iy;

		// update torsion and enforce unit constraint
		xm::vec4 change = .5*dt*xt::linalg::dot(Qmat(),w);
		q += change;
		q *= 1/xt::norm_l2(q)();
	}

	// reset force/torque for next iter
	F.fill(0);
	T.fill(0);
}