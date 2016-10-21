#include <cmath>

#define N 79

constexpr int n_nl = 60;
constexpr int n_timesteps = 10;
constexpr double dt = 0.04;
constexpr double Lx = 79.e3;
constexpr double tol = 1.e-7;
constexpr double adjtol = 1.e-7;
constexpr double ep_glen = 1.e-7;
constexpr double eps = 1.e-5;
constexpr double Aglen = 5.0002e-17;
constexpr double nglen = 3.0;
constexpr double g = 9.81;
constexpr double rhoi = 910.0;
constexpr double rhow = 1035.0;
constexpr double R_bed = -900;
constexpr double beta_const = 5.0;
constexpr double h_left = 1050.0;
constexpr double h_right = 1050.0;
constexpr double dx = Lx / N;
constexpr double PI = 3.14159265358979323844;
scalar r[N];
scalar p[N];
scalar ax[N];

template <typename T>
void stream_vel_init(T*, T*);
template <typename T>
void stream_vel_timedep_stage(T*, T*, T*);
template <typename T>
void stream_vel(T*, T*, T*);
template <typename T>
void stream_vel_taud(T*, T*, T&);
template <typename T>
void stream_vel_visc(T*, T*, T*);
template <typename T>
void stream_vel_assemble(T*, T*, T[N][3]);
template <typename T>
void forward_step(T*, T*, T*);
template <typename T>
void phistage(T*, T*, T*, T*, T*);
template <typename T>
void phi(T*, T*, T*, T*, T*);
template <typename T>
void solve(T*, T*, T[N][3]);
template <typename T>
void solve_lu(T*, T*, T[N][3]);

// here reference on fc mimics the intent(inout) for non-pointer types
template <typename T>
void stream_vel_timedep(T* h, T* u, T& fc) {
  T beta_fric[N];

  for (int i = 0; i < N; i++) {beta_fric[i] = beta_const;}

  stream_vel<T>(u, h, beta_fric);

  stream_vel_timedep_stage<T>(h, u, beta_fric);

  fc = 0; 
  for (int i = 1; i <= N; i++) {
    fc += u[i]*u[i];
  } 

}

template <typename T>
void forward_step(T* h, T* u, T* beta_fric) {
  h[0] = h[0] - dt/dx*u[1]*h[1]; 
  for (int i = 1; i < N; i++) {
    h[i] = h[i] - dt/dx*(u[i+1]*h[i] - u[i]*h[i-1]);
  }
  stream_vel<T>(u, h, beta_fric);
}

template <typename T>
void stream_vel_timedep_stage(T* h, T* u, T* beta_fric) {
  for (int i =0 ; i < n_timesteps; i++ ){
    forward_step<T>(h, u, beta_fric);
  }
}

/*
template <typename T>
void stream_vel_init(T* h, T* beta) {
  for (int i = 0; i < N; i++) {
    beta[i] = beta_const;
    h[i] = h_left + (h_right - h_left) / Lx * ((i+1) - 0.5) * dx;
  }
}
*/

template <typename T>
void stream_vel(T* u, T* h, T* beta_fric) {
  T f[N];
  T utmp[N];
  T b[N];
  T unew[N+1];
  T fend;
  
  // call stream_vel_taud(h, f, fend);
  stream_vel_taud<T>(h, f, fend);

  for (int i = 0; i <= N; i++) {u[i] = 0;}
// ----- driving stress -----

  for (int i = 0; i < N; i++) {
    b[i] = -dx * f[i];
    if (i < N-1) {
      b[i] = b[i] - f[i+1] * dx;
    }
  }
  b[N-1] = b[N-1] + fend;

// --------------------------
  // call phistage(u, unew, b, h, beta_fric, ininloop0)
  phistage<T>(u, unew, b, h, beta_fric); 
  // u = unew;
  for (int i = 0; i <=N; i++) {u[i] = unew[i];}

  // I dont' underwhat is isinloop does in original fortran code
  // But it seems that only one call of the phistage is enough?
}

template <typename T>
void phistage(T* u, T* u_ip1, T* b, T* h, T* beta_fric) {
  phi<T>(u, u_ip1, b, h, beta_fric);
}

template <typename T>
void phi(T* u_i, T* u_ip1, T* b, T* h, T* beta_fric) {
  T nu[N];
  T utmp[N];
  T A[N][3];
  // stream_vel_visc(h, u_i, nu)
  stream_vel_visc<T>(h, u_i, nu);
  
  stream_vel_assemble<T>(nu, beta_fric, A);
  for (int i = 0; i < N; i++) {utmp[i] = 0.0;}
  
  //solve<T>(utmp, b, A);
  solve_lu<T>(utmp, b, A);
  for (int i = 0; i < N; i++) {
    u_ip1[i+1] = utmp[i];
  }
}

template <typename T>
void stream_vel_visc(T* h, T* u, T* nu) {
  T ux, tmp;
  for (int i = 0; i < N; i++) {
    ux = (u[i+1] - u[i]) / dx;
    tmp = ux * ux + ep_glen * ep_glen;
    nu[i] = 0.5*h[i]*pow(Aglen, -1.0/nglen)*pow(tmp, ((1-nglen)/2./nglen));
  }

}

template <typename T>
void stream_vel_assemble(T* nu, T* beta_fric, T A[N][3]) {
  for (int i = 0; i < N; i++) {
    A[i][0] = A[i][2] = 0.0;
    A[i][1] = 4*nu[i]/dx + dx/3.0 * beta_fric[i] * beta_fric[i];
    if (i > 0) {
      A[i][0] = -4*nu[i]/dx + dx/6.0 * beta_fric[i] * beta_fric[i];
    }
    if (i < N-1) {
      A[i][1] = A[i][1] + 4*nu[i+1]/dx+dx/3.0*beta_fric[i+1]*beta_fric[i+1];
      A[i][2] = -4*nu[i+1]/dx + dx/6.0*beta_fric[i+1]*beta_fric[i+1];
    }
  }
}

template <typename T>
void stream_vel_taud(T* h, T* f, T& fend) {
  for (int i = 0; i < N; i++) {
    if (i > 0 && i < N-1) {
      f[i] = rhoi * g * h[i] * (h[i+1] - h[i-1]) / 2.0 / dx;
    } else if (i == 0) {
      f[i] = rhoi * g * h[i] * (h[i+1] - h[i]) / dx;
    } else if (i == N-1) {
      f[i] = rhoi * g * h[i] * (h[i] - h[i-1]) / dx;
    }
  }
  fend = 0.5 * (rhoi * g * h[N-1] * h[N-1] - rhow * g * R_bed*R_bed);
}


template <typename T>
void solve(T* x, T* b, T A[N][3]) {
  T alpha, beta, dp1, dp2, res_init, resid;
  int k_iter = 0;
  for (int i = 0; i < N; i++) {
    r[i] = p[i] = x[i] = 0.0;
  }
  for (int i = 0; i < N; i++) {
    r[i] = b[i] - A[i][1] * x[i];
    if (i > 0) {
      r[i] = r[i] - A[i][0] * x[i-1];
    }
    if (i < N-1) {
      r[i] = r[i] - A[i][2] * x[i+1];
    }
  }  
  dp1 = 0;
  for (int i = 0; i < N; i++) {
    dp1 += r[i]*r[i];
  }
  res_init = sqrt(dp1);
  resid = res_init;
  for (int i = 0; i < N; i++) {
    p[i] = r[i];
  }

  while (k_iter < 200 * N && resid > 1.0e-10*res_init) {
    k_iter = k_iter +1;
    for (int i = 0; i < N; i++) {
      ax[i] = A[i][1] * p[i];
      if (i > 0) {
        ax[i] = ax[i] + A[i][0] * p[i-1];
      }
      if (i < N-1) {
        ax[i] = ax[i] + A[i][2] * p[i+1];
      }
    }
    dp2 = 0;
    for (int i = 0; i < N; i++) {
      dp2 += p[i] * ax[i];
    }
    alpha = dp1 / dp2;

    for (int i = 0; i < N; i++) {
      x[i] = x[i] + alpha*p[i];
      r[i] = r[i] - alpha*ax[i];
    }
    dp2 = dp1; // norm2(r_old);
    dp1 = 0;
    for (int i = 0; i < N; i++) {
      dp1 += r[i] * r[i];
    }
    beta = dp1 / dp2;
    resid = sqrt(dp1);    

    for (int i = 0; i < N; i++) {
      p[i] = r[i] + beta * p[i];
    }
    
  }
}

template <typename T>
void solve_lu(T* x, T* b, T A[N][3]) {
  T d[N];
  T l[N];
  d[0] = sqrt(A[0][1]); l[0] = A[0][2] / d[0];
  for (int i = 1; i < N; i++) {
    d[i] = sqrt(A[i][1] - A[i][0]*A[i][0]/(d[i-1]*d[i-1]));
    l[i] = A[i][2] / d[i];
  }
  // A = L * L^T x = b;
  // L * y = b;
  T y[N];
  for (int i = 0; i < N-1; i++) {
    y[i] = b[i] / d[i];
    b[i+1] -= y[i] * l[i];
  }
  y[N-1] = b[N-1] / d[N-1];
  // L^T x = y;
  for (int i = N-1; i > 0; i--) {
    x[i] = y[i] / d[i];
    y[i-1] -= x[i] * l[i-1];
  }
  x[0] = y[0] / d[0];
}
