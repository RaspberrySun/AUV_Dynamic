import sympy as sp
from sympy import sin, cos, solve
sp.init_printing()

u_ = sp.symbols('u_')
v_ = sp.symbols('v_')
w_ = sp.symbols('w_')
p_ = sp.symbols('p_')
q_ = sp.symbols('q_')
r_ = sp.symbols('r_')
x = sp.symbols('x')
y = sp.symbols('y')
z = sp.symbols('z')
phi = sp.symbols('phi')
theta = sp.symbols('theta')
psi = sp.symbols('psi')
u = sp.symbols('u')
v = sp.symbols('v')
w = sp.symbols('w')
p = sp.symbols('p')
q = sp.symbols('q')
r = sp.symbols('r')
tau_1 = sp.symbols('tau_1')
tau_2 = sp.symbols('tau_2')
tau_3 = sp.symbols('tau_3')
tau_4 = sp.symbols('tau_4')
tau_5 = sp.symbols('tau_5')
tau_6 = sp.symbols('tau_6')

M = sp.Matrix([[187.394, 0, 0, 0, 6.25, 0],
               [0, 187.394, 0, -6.25, 0, 0],
               [0, 0, 187.394, 0, 0, 0],
               [0, -6.25, 0, 4.629, 0, 0],
               [6.25, 0, 0, 0, 4.6290, 0],
               [0, 0, 0, 0, 0, 4.629]])

C_RB = sp.Matrix([[0, 0, 0, 6.25*r, 125*w, -125*v],
                  [0, 0, 0, -125*w, 6.25*r, 125*u],
                  [0, 0, 0, 125*(v-0.05*p), -125*(u+0.05*q), 0],
                  [-6.25*r, 125*w, -125*(v-0.05*p), 0, 4.629*r, -4.629*q],
                  [-125*w, -6.25*r, 125*(u+0.05*q), -4.629*r, 0, 4.629*p],
                  [125*v, -125*u, 0, 4.629*q, -4.629*p, 0]])

C_A = sp.Matrix([[0, 0, 0, 0, 125*w, -125*v],
                 [0, 0, 0, -125*w, 0, 125*u],
                 [0, 0, 0, 125*v, -125*u, 0],
                 [0, 125*w, -125*v, 0, 0, 0],
                 [-125*w, 0, 125*u, 0, 0, 0],
                 [125*v, -125*u, 0, 0, 0, 0]])

D = sp.Matrix([[148*abs(u)+100, 0, 0, 0, 0, 0],
               [0, 148*abs(v)+100, 0, 0, 0, 0],
               [0, 0,  148*abs(w)+100, 0, 0, 0],
               [0, 0, 0, 280*abs(p)+230, 0, 0],
               [0, 0, 0, 0, 280*abs(q)+230, 0],
               [0, 0, 0, 0, 0, 280*abs(r)+230]])

G = sp.Matrix([[2.077*sin(theta)],
               [-2.077*cos(theta)*sin(phi)],
               [-2.077*cos(theta)*cos(phi)],
               [61.3125*cos(theta)*sin(phi)],
               [61.3125*sin(theta)],
               [0]])

V_ = sp.Matrix([[u_],
                [v_],
                [w_],
                [p_],
                [q_],
                [r_]])

V = sp.Matrix([[u], [v], [w], [p], [q], [r]])

Eta = sp.Matrix([[x], [y], [z], [phi], [theta], [psi]])

TAU = sp.Matrix([[tau_1], [tau_2], [tau_3], [tau_4], [tau_5], [tau_6]])

D = M * V_ + (C_RB+C_A)*V + D*V + G - TAU

print(-D)
result = solve([-6.25*p*r - 250*q*w - 6.25*q_ + 250*r*v + tau_1 - u*(148*abs(u) + 100) - 187.394*u_ - 2.077*sin(theta),
           250*p*w + 6.25*p_ - 6.25*q*r - 250*r*u + tau_2 - v*(148*abs(v) + 100) - 187.394*v_ + 2.077*sin(phi)*cos(theta),
           -p*(-6.25*p + 250*v) - q*(-6.25*q - 250*u) + tau_3 - w*(148*abs(w) + 100) - 187.394*w_ + 2.077*cos(phi)*cos(theta),
           -p*(280*abs(p) + 230) - 4.629*p_ + 6.25*r*u + tau_4 - 250*v*w + 6.25*v_ - w*(6.25*p - 250*v) - 61.3125*sin(phi)*cos(theta),
           -q*(280*abs(q) + 230) - 4.629*q_ + 6.25*r*v + tau_5 + 250*u*w - 6.25*u_ - w*(6.25*q + 250*u) - 61.3125*sin(theta),
           -r*(280*abs(r) + 230) - 4.629*r_ + tau_6], [u_, v_, w_, p_, q_, r_])
print(result)
