### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ a0b2fd6d-6ed8-4039-89df-d236e6c6b56b
begin
	using ControlSystems, PlutoUI, Plots, LinearAlgebra
	# Additional packages used in this notebook
end

# ╔═╡ 38173cf4-d94f-11eb-3315-254f36a1f4b9
md"""
# Control Engineering Basics

	Dr. Varodom Toochinda
	Dept. of Mechanical Engineering, Kasetsart University

## Module 5: State Feedback

**Requirement :** Julia with ControlSystems package. To install, in the Julia REPL:

```julia
using Pkg; Pkg.add("ControlSystems")
```

"""

# ╔═╡ 03f53980-7042-4dce-8900-4717013690a7
md"""
This article is contained in Control Engineering Basics study module, which is used as course material for Electrical-Mechanical Manufacturing Engineering (EMME), Department of Mechanical Engineering, Kasetsart University.

### Module Key Study Points
* Understand state-space representation of a system
* How to convert data between state-space and transfer function form
* Design state feedback using simple pole-placement procedure
* Append an integrator to state feedback to eliminate steady-state error
"""

# ╔═╡ 81348b17-d0f4-436f-9713-d347323b89c7
md"""
State feedback control, the topic of this study module, can be thought of as a foundation for the so-called 
*“modern control”* (+) originated since ’60.  In contrast to the frequency domain analysis of the classical 
control theory, modern control theory formulates the time-domain data as a system of first-order differential 
equations in matrix form called *state-space representation*, which is a mathematical model of the given 
system as a set of input, output, and state variable. 

	(+) This terminology is commonly used in the control literature, regardless of referring to a 50-year-old approach as "modern" could sometimes create confusion to a beginner. Indeed, the modern control era is so eternal that  later developments have to be called "post-modern."


First, we give some review on state-space representation, which is essential for the state feedback design 
discussed later on. Let us take our simple DC motor joint model as an example. A dynamic equation that governs the 
motion can be written as

```math 
J\ddot{\theta}(t) + B\dot{\theta}(t) = u(t) \tag{1}
```

To simplify the notation, the time dependent is omitted. Define the system states as the joint position and velocity

```math
\begin{matrix}
x_1 &= \theta  \\
x_2 &= \dot{\theta}
\end{matrix} \tag{2}
```

By using these state variables, (1) can be rewritten as a system of first order differential equations

```math
\begin{eqnarray}
\dot{x}_1 &=& x_2  \\
\dot{x}_2 &=& -\frac{B}{J}x_2 + \frac{1}{J}u(t)
\end{eqnarray} \tag{3}
```
or in matrix form as


```math 
\left[ \begin{array}{c}
\dot{x}_1 \\
\dot{x}_2 
\end{array}  \right] = 
\left[ \begin{array} {cc}
0 & 1 \\
0 & -\frac{B}{J}
\end{array} \right]
\left[ \begin{array}{c}
x_1 \\
x_2
\end{array} \right] + 
\left[ \begin{array}{c}
0 \\
\frac{1}{J} 
\end{array} \right] u \tag{4}
```
and the output equation, with $y = \theta$

```math
y = \left[ \begin{array}{cc}
1 & 0
\end{array} \right] \left[ \begin{array}{c}
x_1 \\
x_2
\end{array} \right] \tag{5}
```
In general, a linear time-invariant (LTI) system can be described as

```math
\dot{x} = Ax + Bu \tag{6}
```

```math
y = Cx + Du \tag{7}
```
where $x,u,y$ represent the state, input, and output vectors, respectively.  (6) is called a 
*state equation*, and (7) an *output equation*. Note that this representation is not unique, 
but depends on how the states are defined. Moreover, though in our robot joint example the states are conveniently 
joint position and velocity, for the general case the states might not have any physical meaning, and hence could 
not be measured. 

For a system represented by (6) and (7), it is easy to show that the corresponding transfer function equals

```math 
P(s) = C(sI - A)^{-1}B + D \tag{8}
```
To convert between state space and transfer function with Python Control library, use ss2tf() and tf2ss(). 
For example, for a plant transfer function 

```math
P(s) = \frac{1}{10s^2 + 0.1s}
```
"""

# ╔═╡ e25e7e4d-ca21-4c2a-bc57-7e0a645053ff
begin
	s = tf("s")
	P = 1/(10s^2+0.1s)
end

# ╔═╡ 9ac2a271-fef7-4a56-847b-5cc87230e94f
md"""
This can be converted to state-space form by
"""

# ╔═╡ e55c42dc-8ed2-4036-a577-70635a0367c3
Pss = ss(P)

# ╔═╡ a86dba5a-50f9-4778-842f-5591c1f6bd24
md"""
Note some slight difference in some elements of $B, C$ from (3),(4). This result reveals that state-space representation is not unique. 

To convert back to transfer function, use tf() function on Pss. We see that some normalization is applied such that the coefficient of highest order of $s$  in the denominator is one.
"""

# ╔═╡ 7441c155-5dc2-4e6e-963e-bb4026fef6d4
P₁ = tf(Pss)

# ╔═╡ df4bff16-f2c1-40ab-a986-5faa6936d84f
md"""
Note that a transfer function for a physical system must be strictly proper; i.e., its frequency response must go to
zero as the frequency approaches infinity. (No system could have unlimited bandwidth in reality.) This implies its state-space 
representation must have zero $D$ matrix. 

### State Feedback Control

Obviously, a state feedback control is feasible in practice when all states are measurable, such as the robot joint 
dynamics in (4) with joint position and velocity as state variables. State feedback design can be performed with a scheme 
known as *pole placement*, or more systematic way, using *Ackerman’s formula*. 

The idea of pole placement is simple. Specify pole locations of the closed-loop system that yields desired stability
and performance criteria, for instance, in some shaded area shown in Figure 1. Then choose the state feedback control gains 
to move the closed-loop poles to such locations. A necessary condition is that the plant must be stabilizable. 
The details can be studied from most undergraduate control textbooks. 

!["Fig1"](https://drive.google.com/uc?id=1OIVp0VUnhaXckVQufic17-aA13A-8DwK)

Figure 1 desired closed-loop pole locations

The state feedback controller is simply a vector of gains $K = \left[k_1,  k_2,  \ldots , k_n \right]^T$ connecting the states  to the plant input. So, for a set of specified closed-loop poles, the design goal is to 
compute $k_i, i = 1, \ldots ,n$. 

Assume for the moment that the command input is zero. As shown in Figure 2, we have at the plant input

!["Fig2"](https://drive.google.com/uc?id=14Mok9cUbcqwNSCNjo1TMrSJLX2jJHsJf)

Figure 2 block diagram of state feedback control

```math
u = -Kx \tag{9}
```

and the closed-loop state equation

```math
\dot{x} = (A - BK)x \tag{10}
```

The closed-loop poles can be computed from

```math
det(sI - A + BK) = 0 \tag{11}
```	
Meanwhile, specifying the closed-loop poles $p_i, i = 1, \ldots , n$ yields the characteristic polynomial

```math
\alpha(s) = (s - p_1)(s - p_2) \ldots (s - p_n) \tag{12}
```
Hence, we can compare (11) and (12) to solve for $K$ manually,  which could be tedious for higher order equations.
Python control library provides a convenient function place() to solve for $K$, given the $A, B$ matrices and a 
vector of desired poles as arguments.

**Example 1:** Let us design a state feedback control for the simple robot joint described by (4), (5) 
with $J = 10, B = 0.1$. Specify the desired properties of closed-loop system as follows:

1. overshoot less than 5%
2. rise time less than 0.1 sec

By standard analysis of 2nd order system, we have that the specification 1 translates to damping ratio $\zeta \ge 0.7$, and using the relation $t_r = 1.8/\omega_n$, we have for specification 2 that $\omega_n \ge 18$ rad/s.  Substituting these two values to the closed-loop characteristic polynomials yields 

```math 
\Lambda(s) = s^2 + 2\zeta\omega_ns + \omega_n^2 = s^2 + 25.2s + 324 \tag{12}
```
with poles at $-12.6 \pm 12.8546i$. The above procedure can be carried out by the following code
"""

# ╔═╡ c0e42850-d9f3-4399-97a2-03c0ee75c103
begin
	ζ = 0.7
	ωₙ = 18
	Λ = s^2 + 2ζ*ωₙ*s + ωₙ^2
	desired_poles = tzero(Λ)
end

# ╔═╡ 4360f113-8715-441f-869f-b678b2fac6ea
md"""
Construct the plant. It can be created first in transfer function form and converted to state-space. For this problem, we
choose to create the $A, B, C, D$ matrix directly using (4),(5).
"""

# ╔═╡ dc5cde72-50e7-4372-8c44-0b114c8d8f3a
begin
	A = [0.0 1.0;0.0 -0.01]
	B = [0.0; 0.1]
	C = [1.0 0.0]
	D = [0]
	P_ss = ss(A,B,C,D)
end

# ╔═╡ ed328876-2f0c-46eb-973a-f31bc00ed4d1
md"""
Before performing state feedback control, it is safe to test whether the plant is controllable with the function ctrb(). The 
controllability matrix should have full rank, in this case 2.

	Note : rank() requires LinearAlgebra package
"""

# ╔═╡ f6e0fd19-459a-47db-96ae-d07bafc1a185
rank(ctrb(P_ss))

# ╔═╡ 6a385b9b-9fac-460d-8d8b-b18a34658a8a
md"""
Now we can use place() to compute the state feedback gain.
"""

# ╔═╡ 3900ab40-5e6a-46dd-ae0b-0c091d3b2974
K = real(place(P_ss, desired_poles))

# ╔═╡ 22e95c3b-e1ea-4d68-a4fe-5bec6860d5a6
md"""
Check the poles of the closed-loop system to verify that they are in the desired locations
"""

# ╔═╡ 7c29e1b9-be44-4581-8dec-36853050f7b9
begin
	T_ss = ss(A-B*K, B, C, D)
	pole(T_ss)
end

# ╔═╡ 8fba9dd8-e5b9-4718-b05f-893d55739351
md"""
The original simulation model was built using Scilab/xcos as shown in Figure 3. Notice in the diagram that the plant 
is conveniently represented in transfer function form since the joint velocity and angle can be accessed. 

!["Fig 3"](https://drive.google.com/uc?id=1ZuLvaxvGSIqlRxfHX0AX3I-JG7PkHcfW)

Figure 3 Scilab/xcos diagram to simulate state feedback control

Also, in computing the state feedback gains, we do not take into consideration the command input. Hence the step 
response will have nonzero steady-state error that needs to be compensated with a feedforward gain. The DC gain can 
be computed using dcgain(). Therefore, the feedforward gain is 

"""

# ╔═╡ cd88d8dd-2a59-4821-8494-115735e4daa1
ffgain = real((1 ./dcgain(T_ss))[1])

# ╔═╡ 32da649d-1e6e-4fa3-a8a2-dce45ffc0090
md"""
With the step disturbance of magnitude 0.1 injecting to the system at time t = 1 sec, the simulation yields the 
step response in Figure 4. We see that the transient period conforms to the desired spec; i.e., 1. overshoot less 
than 5%. 2. rise time less than 0.1 sec. However, the closed-loop system cannot get rid of the constant disturbance 
after t = 1 sec. This result is predictable, because the state feedback is just a pair of static gains with no 
dynamics to compensate the disturbance entering at the plant input.

!["Fig 4"](https://drive.google.com/uc?id=1pQ_GNQ6ElegD3TAG1z02cqpeog8u05IP)

Figure 4 step response from the model in Figure 3

To perform this simulation using Julia with ControlSystems package, we simply create an input signal that is a sum of 
the reference command and the input disturbance at t = 1 sec. This somehow gives us more insight to the problem. The 
closed-loop system can never get rid of the disturbance because it is seen by the feedback system as the reference 
command changing its level.
"""

# ╔═╡ 3a956dd5-7aa1-496a-bada-541a1ef1091c
begin
	tvec = collect(0:0.001:2)
	midpts = trunc(Int,length(tvec)/2)
	dist_level = 0.1
	r = ones(size(tvec))
	rc = ffgain*r
	d = dist_level*ffgain*r
	d[1:midpts] .= 0.0
	u = rc - d
	y, t, x = lsim(T_ss, u, tvec, method=:zoh)
	plot(t,r,label = "step ref")
	plot!(t,y, label= "plant output",xlabel="time (sec)",ylabel="y(t)",title="Step response of state feedback control",legend=:bottomright)
	
	
end

# ╔═╡ c5dda403-25f5-4e88-903d-a33bf53af762
md"""
**Example 2:** From the PID example discussed earlier, we show the advantage of integral term in eliminating 
the steady-state error. This principle can be applied to the state feedback scheme by augmenting an integrator as 
shown by the simulation diagram in Figure 5. There exist some systematic procedure to augment an integrator and design
all gains simultaneously, but that means the second-order relationship in the previous example is no longer 
applicable. So in this example we still use the pole-placement for second-order system as before, and then adjust 
the integral gain afterwards to achieve the desired response.

!["Fig 5"](https://drive.google.com/uc?id=1wIUhJGgMMsTe75qYS76ufcDUxTf99-m-)

Figure 5 Xcos simulation model for the state feedback with integrator control

Notice that the integrator replaces the DC gain compensation in Figure 1 to correct the response to the target 
steady-state value. Using the same state feedback gains results in slower transient response than Figure 4, 
so we redesign the pole-placement step by increasing $\omega_n$   to 40 rad/s.

**Note :** to be more flexible, $\omega_n$ is bound to a slider with default value of 40 rad/s
"""

# ╔═╡ 4950cce7-43d9-4469-9d5c-39b60a9b3a12
md"""
This yields a new pair of state feedback gains.
"""

# ╔═╡ 0baa19ce-5450-442b-abfe-8e449cc46291
md"""
and the feedforward gain
"""

# ╔═╡ 54007e31-3f5b-4e87-8e98-c02a065f146a
md"""
is applied as compensation to the disturbance input to yield the same level  as in previous example; 
i.e., $d = 0.1$. By experimenting with the integral gain, we select $K_i = 200000$ , which gives the response as 
in Figure 6. The rise time and overshoot satisfy the specifications, while the system recovers to the desired value 
after the disturbance is applied at t = 1 sec. 

!["Fig 6"](https://drive.google.com/uc?id=1G7I8jsn3zJL_Gnl6PzH7RLoIjexv6QYA)

Figure 6 responses from the Xcos model in Figure 5

To simulate this using Julia with ControlSysems package, the integrator have to be augmented to the state-space 
representation. From a general state feedback with integrator diagram in Figure 7, it is straightforward to show 
that the open-loop system is described by

```math 
\left[ \begin{array}{c}
\dot{x} \\
\dot{w} 
\end{array}  \right] = 
\left[ \begin{array} {cc}
A & 0 \\
c & 0
\end{array} \right]
\left[ \begin{array}{c}
x \\
w
\end{array} \right] + 
\left[ \begin{array}{c}
b \\
0 
\end{array} \right] (u + d) + 
\left[ \begin{array}{c}
0 \\
-1 
\end{array} \right] r \tag{13}
```

!["Fig 7"](https://drive.google.com/uc?id=1YzlfVXnCv_luVi8uOagtp779hQA4sekS)

Figure 7 general block diagram of state feedback with integrator control

Applying the control law

```math
u = -Kx -k_iw \tag{14}
```

results in the closed-loop system

```math
\left[ \begin{array}{c}
\dot{x} \\
\dot{w} 
\end{array}  \right] = 
\left[ \begin{array} {cc}
A-bK & -bk_i \\
c & 0
\end{array} \right]
\left[ \begin{array}{c}
x \\
w
\end{array} \right] + 
\left[ \begin{array}{c}
b \\
0 
\end{array} \right] d + 
\left[ \begin{array}{c}
0 \\
-1 
\end{array} \right] r \tag{15}
```

Create this system and do the simulation.
"""

# ╔═╡ 90552203-a914-4e91-b436-0c3ccf5f255d
md"""
Kᵢ = $(@bind Kᵢ Slider(1:1000000, show_value=true, default=200000))
, ωₙ₁ = $(@bind ωₙ₁ Slider(10:100, show_value=true, default=40)) rad/s

"""

# ╔═╡ d07b6d1e-7184-465e-8926-2d37e7176470
begin
	# ζ = 0.7
	#ωₙ₁ = 40  # use Slider below to adjust 
	Λ₁ = s^2 + 2ζ*ωₙ₁*s + ωₙ₁^2
	desired_poles_1 = tzero(Λ₁)
end

# ╔═╡ ec3b3cfc-8d9d-460e-aa51-7a7d4d765f6d
K₁ = real(place(P_ss, desired_poles_1))

# ╔═╡ b650b8b5-0dc5-4e2e-b36d-201b806aac42
begin
	Tss_1 = ss(A-B*K₁, B, C, D)
	ffgain_1 = real((1 ./dcgain(Tss_1))[1])
end

# ╔═╡ cd197e58-f07b-4417-80e8-7e612f48dee4
begin
	A₁ = [A - B*K₁ -Kᵢ*B]
	Cₐ = [C [0]]
	Aₐ = [A₁;Cₐ]
	B₁ = [B; 0]
	B₂ = [0;0;-1]
	Bₐ = [B₁ B₂]
	Dₐ = zeros(1,2)
	Tssₐ = ss(Aₐ,Bₐ,Cₐ,Dₐ)
end

# ╔═╡ 93d7de76-3170-4be6-a99a-3b82d6fbf9ab
begin
	#tvec = collect(0:0.001:2)
	#midpts = trunc(Int,length(tvec)/2)
	#dist_level = 0.1
	#r = ones(size(tvec))
	d₁ = -dist_level*ffgain_1*r
	d₁[1:midpts] .= 0.0
	u₁ = [d₁ r]
	yₐ, tₐ, xₐ = lsim(Tssₐ, u₁, tvec, method=:zoh)
	plot(tₐ,r,label = "step ref")
	plot!(tₐ,yₐ, label= "plant output",xlabel="time (sec)",ylabel="y(t)",title="Step response of state feedback control",legend=:bottomright)
	
	
end

# ╔═╡ 4b5c72ea-313e-4b79-971b-7eb67a62889b
md"""
### Summary

In this module, we discuss state-space representation, using the robot joint driven by DC motor model as an 
example. The two state variables are the joint position and velocity, which are measurable in a real application. 
Hence the state feedback design scheme is suitable for this system. Joint position is normally obtained from an 
encoder using hardware or software readouts. Joint velocity may be measured via a tachometer, or obtained indirectly 
by counting encoder pulses per known time period. Finally, we show how to append an integrator to a state feedback 
design to eliminate steady-state error.  

### References

1. V.Toochinda. ["Robot Analysis and Control with Scilab and RTSX"](http://dewninja.blogspot.com/p/robot-analysis-and-control-with-scilab.html) . Mushin Dynamics, 2014. 
"""

# ╔═╡ a2c4f79b-2961-4a6e-9ec9-25283cc95ca8


# ╔═╡ efb833a4-17fc-4569-b4e0-b9cba88f53cd
md"""
**Last updated :** July 22, 2021
"""

# ╔═╡ e988f4b8-3c94-4916-ba0c-5d1c8d3d6430
md"""
!["dewninja"](https://drive.google.com/thumbnail?id=13bzT7Rmy3bzvE7TiS0yfQo94kpxMuipF)

	dew.ninja
	Copyright 2021
"""

# ╔═╡ Cell order:
# ╟─38173cf4-d94f-11eb-3315-254f36a1f4b9
# ╟─03f53980-7042-4dce-8900-4717013690a7
# ╠═a0b2fd6d-6ed8-4039-89df-d236e6c6b56b
# ╟─81348b17-d0f4-436f-9713-d347323b89c7
# ╠═e25e7e4d-ca21-4c2a-bc57-7e0a645053ff
# ╟─9ac2a271-fef7-4a56-847b-5cc87230e94f
# ╠═e55c42dc-8ed2-4036-a577-70635a0367c3
# ╟─a86dba5a-50f9-4778-842f-5591c1f6bd24
# ╠═7441c155-5dc2-4e6e-963e-bb4026fef6d4
# ╟─df4bff16-f2c1-40ab-a986-5faa6936d84f
# ╠═c0e42850-d9f3-4399-97a2-03c0ee75c103
# ╟─4360f113-8715-441f-869f-b678b2fac6ea
# ╠═dc5cde72-50e7-4372-8c44-0b114c8d8f3a
# ╟─ed328876-2f0c-46eb-973a-f31bc00ed4d1
# ╠═f6e0fd19-459a-47db-96ae-d07bafc1a185
# ╟─6a385b9b-9fac-460d-8d8b-b18a34658a8a
# ╠═3900ab40-5e6a-46dd-ae0b-0c091d3b2974
# ╟─22e95c3b-e1ea-4d68-a4fe-5bec6860d5a6
# ╠═7c29e1b9-be44-4581-8dec-36853050f7b9
# ╟─8fba9dd8-e5b9-4718-b05f-893d55739351
# ╠═cd88d8dd-2a59-4821-8494-115735e4daa1
# ╟─32da649d-1e6e-4fa3-a8a2-dce45ffc0090
# ╠═3a956dd5-7aa1-496a-bada-541a1ef1091c
# ╟─c5dda403-25f5-4e88-903d-a33bf53af762
# ╠═d07b6d1e-7184-465e-8926-2d37e7176470
# ╟─4950cce7-43d9-4469-9d5c-39b60a9b3a12
# ╠═ec3b3cfc-8d9d-460e-aa51-7a7d4d765f6d
# ╟─0baa19ce-5450-442b-abfe-8e449cc46291
# ╠═b650b8b5-0dc5-4e2e-b36d-201b806aac42
# ╟─54007e31-3f5b-4e87-8e98-c02a065f146a
# ╠═cd197e58-f07b-4417-80e8-7e612f48dee4
# ╟─90552203-a914-4e91-b436-0c3ccf5f255d
# ╠═93d7de76-3170-4be6-a99a-3b82d6fbf9ab
# ╟─4b5c72ea-313e-4b79-971b-7eb67a62889b
# ╠═a2c4f79b-2961-4a6e-9ec9-25283cc95ca8
# ╟─efb833a4-17fc-4569-b4e0-b9cba88f53cd
# ╟─e988f4b8-3c94-4916-ba0c-5d1c8d3d6430
