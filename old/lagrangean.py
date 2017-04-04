import numpy as np

from kpop.admixture.mixture_algorithm import AdmixtureCalculation


class LagrangeanMinimizer(AdmixtureCalculation):

    def __init__(self, data, k, freqs=None, q=None, gamma_q=0.5, gamma_f=0.5, mass_f=1, mass_q=100, baumgarte_gamma=0.1, vq0_mul=0.00, vf0_mul=0.0):
        super().__init__(data, k, freqs, q)
        self.vel_f = -np.array([self.jac_fj(j, self.freqs[j])
                                for j in range(self.num_loci)]) * vf0_mul
        self.vel_q = -np.array([self.jac_qi(i, self.q[i])
                                for i in range(self.num_ind)]) * vq0_mul
        self.vel_baumgarte_q = self.q * 0.0
        self.gamma_q = gamma_q
        self.gamma_f = gamma_f
        self.mass_f = np.zeros(self.num_loci, dtype=float) + mass_f
        self.mass_q = np.zeros(self.num_ind, dtype=float) + mass_q
        self.baumgarte_gamma = baumgarte_gamma
        self.vels_q = []
        self.vels_f = []
        self.forces_q = []
        self.forces_f = []
        self.vels_diff_q = []
        self.vels_diff_f = []

    def advance(self, dt):
        freqs = self.freqs
        q = self.q
        vel_f = self.vel_f
        vel_q = self.vel_q

        # Compute forces
        force_f = -np.array([self.jac_fj(j, freqs[j])
                             for j in range(self.num_loci)])
        force_q = -np.array([self.jac_qi(i, q[i])
                             for i in range(self.num_ind)])
        vf = np.sqrt((m.vel_f**2).sum(1))
        vq = np.sqrt((m.vel_q**2).sum(1))

        force_f -= self.mass_f[:, None] * self.gamma_f * vel_f
        force_q -= self.mass_q[:, None] * self.gamma_q * vel_q

        vel_f_ = vel_f + force_f * dt / self.mass_f[:, None]
        vel_q_ = vel_q + force_q * dt / self.mass_q[:, None]

        impulse_f = force_f * dt / self.mass_f[:, None]
        impulse_q = force_q * dt / self.mass_q[:, None]

        # Apply constraints in f
        for j, fj in enumerate(freqs):
            for k, f_jk in enumerate(fj):
                if (f_jk <= 0 and vel_f_[j, k] < 0) or (f_jk >= 1 and vel_f_[j, k] > 0):
                    impulse_f[j, k] -= vel_f_[j, k]

        # Apply constraints in q
        for i, qi in enumerate(q):
            J = np.ones(self.k, dtype=float)
            for k, q_ik in enumerate(qi):
                if (q_ik <= 0 and vel_q_[i, k] <= 0) or (q_ik >= 1 and vel_q_[i, k] >= 0):
                    impulse_q -= vel_q_[i, k]
                    J[k] = 0
            if J.sum() != 3:
                print(i, qi)
            lambd = -np.dot(J, qi)
            impulse_q += dt * lambd * J

        # Evolve positions
        vel_f += impulse_f
        vel_q += impulse_q
        freqs += vel_f * dt
        q += vel_q * dt
        vf_ = np.sqrt((m.vel_f**2).sum(1)).mean()
        vq_ = np.sqrt((m.vel_q**2).sum(1)).mean()
        f_f = abs(impulse_f).mean() / dt
        f_q = abs(impulse_q).mean() / dt
        print('vels-post', vf_, vq_)
        print('forces', f_f, f_q)
        self.forces_f.append(f_f)
        self.forces_q.append(f_q)
        self.vels_diff_f.append(vf.mean() - vf_)
        self.vels_diff_q.append(vq.mean() - vq_)
        self.vels_f.append(vq)
        self.vels_q.append(vq)

        # Enforce constraints
        for j, fj in enumerate(freqs):
            for k, f_jk in enumerate(fj):
                if f_jk < 0:
                    fj[k] = 0
                elif f_jk > 1:
                    fj[k] = 1

        for i, qi in enumerate(q):
            for k, q_ik in enumerate(qi):
                if q_ik < 0:
                    qi[k] = 0
                elif q_ik > 1:
                    qi[k] = 1
            qi /= qi.sum()


def pos_update(dt, check=True, baumgarte=False, normalize=False):
    global F, Q, VQ, VF

    F += VF * dt
    Q += VQ * dt

    if baumgarte:
        error = (Q.sum(1) - 1) / 3
        Q -= 0.9 * error[:, None]  # baumgarte correction

    if check:
        VQ = np.where((Q < 0) * (VQ < 0) + (Q > 1) * (VQ > 0), -VQ, VQ)
        VF = np.where((F < 0) * (VF < 0) + (F > 1) * (VF > 0), -VF, VF)
        Q = np.where(Q < 0, 0, Q)
        Q = np.where(Q > 1, 1, Q)
        F = np.where(F < 0, 0, F)
        F = np.where(F > 1, 1, F)

    if normalize:
        Q /= Q.sum(1)[:, None]


def pos_update_r(dt, check=True):
    global F, Q, VQ, VF

    F += VF * dt
    Qr += VQr * dt
    Q[:, :2] = Qr
    Q[:, 2] = Qr.sum(1)

    if check:
        VQr = np.where((Qr < 0) * (VQr < 0) + (Qr > 1) * (VQr > 0), -VQr, VQr)
        VFr = np.where((F < 0) * (VF < 0) + (F > 1) * (VF > 0), -VF, VF)
        Qr = np.where(Qr < 0, 0, Qr)
        Qr = np.where(Qr > 1, 1, Qr)
        F = np.where(F < 0, 0, F)
        F = np.where(F > 1, 1, F)

    if normalize:
        Q /= Q.sum(1)[:, None]


def vel_update(dt, check=True, normalize=False):
    global VF, VQ, impulseF, impulseQ

    Fq, Ff = forces(Q, F)
    impulseF = Ff * dt
    impulseQ = Fq * dt

    if normalize:
        Jn = np.array([1., 1., 1.])
        lambd = -np.dot(VQ + impulseQ, Jn) * (mass_q / 3)
        impulseQ = lambd[:, None] * Jn[None, :]

    VF += impulseF / mass_f
    VQ += impulseQ / mass_q
    Jn = np.ones(K)
    err = Q.sum(1) - 1

    if check:
        VF[:] = np.where((F <= 0) * (F >= 1), 0, VF)
        VQ[:] = np.where((Q <= 0) * (Q >= 1), 0, VQ)


def simulation():
    Q, F = QQ.copy(), FF.copy()
    Q = Q * 0 + 1 / 3
    F = F * 0 + 0.5
    VQ, VF = initial_direction(1.5)
    dt = (1e-4)
    mass_q, mass_f = initial_mass()
    mass_q = mass_f = 1
    times = []
    t, tmax = -dt, 0.22
    reset_tracker()
    while t < tmax:
        t += dt

        vel_update(dt, check=True, normalize=True)
        pos_update(dt, check=True, baumgarte=False, normalize=False)

        track()
        times.append(t)
        if energy_v and np.isnan(energy_v[-1]):
            break


Fcomp = FF.copy()
FBAR_CACHE = G + 0.0
FBAR_CACHE_ = G + 0.0
FF_CACHE = FF * 0.0
FQ1_CACHE = QQ * 0.0
FQ2_CACHE = QQ * 0.0


def forces_python(Q, F):
    global fbar
    F_comp = 1 - F
    fbar = np.dot(Q, F.T)
    fbar_comp = np.dot(Q, F_comp.T)
    chi1 = G / (fbar + e)
    chi2 = (2 - G) / (fbar_comp + e)
    Ff = np.dot((chi1 - chi2).T, Q)
    Fq = np.dot(chi1, F) + np.dot(chi2, F_comp)
    return Fq, Ff


def forces(Q, F):
    return _forces(Q, F, G)
    return _force_q(Q, F, G), _force_f(Q, F, G)


if __name__ == '__main__':
    from kpop import Population
    N = 10
    J = 100
    pA = Population.make_random(N, J, label='A', seed=1)
    pB = Population.make_random(N, J, label='B', seed=2)
    pC = Population.make_random(N, J, label='C', seed=3)
    pop = pA + pB + pC
    m = LagrangeanMinimizer(pop, 3, vf0_mul=0, vq0_mul=0, gamma_q=4, gamma_f=2)

    for _ in range(10):
        m.optimize_em_acc()
        f0 = m.log_likelihood()
        print('EM: log-like: %s' % f0)

    ll = []
    dll = []
    eK = []
    eT = []
    for _ in range(500):
        print('-' * 80)
        L0 = m.log_likelihood()
        m.advance(5e-4)
        L1 = m.log_likelihood()
        print('log-like: %s (delta=%s)' % (L1, L1 - L0))
        ll.append(L1)
        dll.append(L1 - L0)
        eK.append(((m.vel_f**2).sum(1) * m.mass_f / 2).sum() +
                  ((m.vel_q**2).sum(1) * m.mass_q / 2).sum())
        eT.append(eK[-1] + L1)

    from matplotlib import pyplot as plt
    #plt.plot(m.vels_diff_f); plt.title('vels diff (f)'); plt.admix()
    #plt.plot(m.vels_diff_q); plt.title('vels diff (q)'); plt.admix()
    plt.plot(m.vels_f)
    plt.title('mean vels (f)')
    plt.show()
    plt.plot(m.vels_q)
    plt.title('mean vels (q)')
    plt.show()
    plt.plot(m.forces_f)
    plt.title('forces (f)')
    plt.show()
    plt.plot(m.forces_q)
    plt.title('forces (q)')
    plt.show()
    plt.plot(ll)
    plt.title('log-like')
    plt.show()
    plt.plot(dll)
    plt.title('log-like delta')
    plt.show()
    plt.plot(eK)
    plt.title('kinetic energy')
    plt.show()
    plt.plot(eT)
    plt.title('total energy')
    plt.show()
