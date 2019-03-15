from scipy.integrate import odeint
import numpy as np
from math import cos, sin, tan
import matplotlib.pyplot as plt


def six_dof(states, t, tau_1, tau_2, tau_3, tau_4, tau_5, tau_6):
    # 计算12个状态变量x, y, z, phi, theta, psi, u, v, w, p, q, r
    x, y, z, phi, theta, psi, u, v, w, p, q, r = states
    return np.array([cos(psi)*cos(theta)*u+(cos(psi)*sin(theta)*sin(phi))*v + (cos(psi)*sin(theta)*cos(phi)+sin(phi)*sin(psi))*w,
                     sin(psi)*cos(theta)*u+(sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi))*v+(sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi))*w,
                     -sin(theta)*u+cos(theta)*sin(phi)*v+cos(theta)*cos(phi)*w,
                     p + tan(theta)*sin(phi)*q + cos(phi)*tan(theta)*r,
                     cos(phi)*q-sin(phi)*r,
                     sin(phi)/cos(theta)*q - cos(phi)/cos(theta)*r,
                     -0.0349249123769635 * p * r - 1.34984145028355 * q * w + 2.11254600681568 * q * abs(q) + 1.73530564845574 * q + 1.34984145028355 * r * v + 0.00558798598031417 * tau_1 - 0.00754480716719886 * tau_5 - 0.827021925086497 * u * abs(u) - 0.558798598031417 * u + 0.450984742557768 * sin(theta),
                     1.34984145028355 * p * w - 2.11254600681568 * p * abs(p) - 1.73530564845574 * p - 0.0349249123769635 * q * r - 1.34984145028355 * r * u + 0.00558798598031417 * tau_2 + 0.00754480716719886 * tau_4 - 0.827021925086497 * v * abs(v) - 0.558798598031417 * v - 0.450984742557768 * sin(phi) * cos(theta),
                     0.0333521884371965 * p ** 2 - 1.33408753748786 * p * v + 0.0333521884371965 * q ** 2 + 1.33408753748786 * q * u + 0.00533635014995144 * tau_3 - 0.789779822192813 * w * abs(w) - 0.533635014995144 * w + 0.0110835992614491 * cos(phi) * cos(theta),
                     0.472350197509652 * p * w - 63.3405514241949 * p * abs(p) - 52.0297386698743 * p - 0.0471550447949929 * q * r - 0.472350197509652 * r * u + 0.00754480716719886 * tau_2 + 0.22621625508641 * tau_4 - 1.11663146074543 * v * abs(v) - 0.754480716719886 * v - 13.8542135754993 * sin(phi) * cos(theta),
                     0.0471550447949929 * p * r + 0.472350197509652 * q * w - 63.3405514241949 * q * abs(q) - 52.0297386698743 * q - 0.472350197509652 * r * v - 0.00754480716719886 * tau_1 + 0.22621625508641 * tau_5 + 1.11663146074543 * u * abs(u) + 0.754480716719886 * u - 13.8542135754993 * sin(theta),
                     -60.4882263987902 * r * abs(r) - 49.6867573990063 * r + 0.216029379995679 * tau_6])


t = np.arange(0, 30, 0.01)

track = odeint(six_dof, (0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00), t, args=(0, 0, 100, 0, 0, 0))

x = np.arange(0, 30, 0.01)
plt.plot(x, track[:, -4])
plt.xlabel('Time')
plt.ylabel('ZSpeed')
plt.title("ZForce-100")
plt.show()
