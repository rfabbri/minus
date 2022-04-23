#ifndef internal_util_h_
#define internal_util_h_

#include <random>

namespace MiNuS {
  
template <typename F>
using C = typename std::complex<F>;

template <unsigned N, typename F>
struct minus_array { // Speed critical -----------------------------------------
  static inline void 
  multiply_scalar_to_self(C<F> *__restrict a, C<F> b)
  {
    for (unsigned i = 0; i < N; ++i, ++a) *a = *a * b;
  }

  static inline void
  negate_self(C<F> * __restrict a)
  {
    for (unsigned i = 0; i < N; ++i, ++a) *a = -*a;
  }

  static inline void 
  multiply_self(C<F> * __restrict a, const C<F> * __restrict b)
  {
    for (unsigned int i=0; i < N; ++i,++a,++b) *a *= *b;
  }

  static inline void 
  add_to_self(C<F> * __restrict a, const C<F> * __restrict b)
  {
    for (unsigned int i=0; i < N; ++i,++a,++b) *a += *b;
  }

  static inline void 
  add_scalar_to_self(C<F> * __restrict a, C<F> b)
  {
    for (unsigned int i=0; i < N; ++i,++a) *a += b;
  }

  static inline void 
  copy(const C<F> * __restrict a, C<F> * __restrict b)
  {
    memcpy(b, a, N*sizeof(C<F>));
  }

  static inline F
  norm2(const C<F> *__restrict a)
  {
    F val = 0;
    C<F> const* __restrict end = a+N;
    while (a != end) val += std::norm(*a++);
    return val;
  }
  
  // Get the real part of the solution vector s (e.g., member x of struct
  // solution).
  //  
  // rs: real solution; holds solution R12, t12, R13, T13 row-major
  //
  // \returns true if the solution is nearly real, false otherwise
  //
  // Not speed critical.
  static inline bool
  get_real(const C<F> s[N], F rs[N])
  {
    // Hongyi function realSolutions = parseSolutionString(output)
    // solutions = reshape(solutions,[14,length(solutions)/14]);
    static constexpr double eps = 10e-3;
    
    // TODO[improvement]
    // Fancy way to coNrt to real is to check if the complex number is close to
    // horizontal then get absolute value.
    /*
    for (unsigned var = 0; var < N; ++var)  // differs from Hongyi criterion
      if (s->x[var].real() < eps && s->x[var].real() >= eps
          || std::abs(std::tan(std::arg(s->x[var].imag()))) >= eps)
        return false;
    
    F real_solution[N];
    for (unsigned var = 0; var < N; ++var) 
      real_solution[var] = ((s->x[var].real() >= 0) ? 1 : -1) * std::abs(s->x[var]);
    */
    unsigned var = 0;
    for (; var < N; ++var)
      if (std::abs(s[var].imag()) >= eps) return false;
    for (var = N-1; var != (unsigned)-1; --var) 
      rs[var] = s[var].real();

    // quat12 rs(0:3), quat12 rs(4:7)
    //  T12 = solutions(9:11);
    //  T13 = solutions(12:14);
    //  R12 = quat2rotm(transpose(quat12));
    //  R13 = quat2rotm(transpose(quat13));
    return true;
  }
};

// Functions over 3 dimensions
template <typename F>
struct minus_3d {
  static inline void
  cross(const C<F> v1[3], const C<F> v2[3], C<F> r[3])
  {
    r[0] = v1[1] * v2[2] - v1[2] * v2[1];
    r[1] = v1[2] * v2[0] - v1[0] * v2[2];
    r[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }
  
  static inline C<F>
  dot(const C<F> v1[3], const C<F> v2[3]) { return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
  
  static inline void
  cross(const F v1[3], const F v2[3], F r[3])
  {
    r[0] = v1[1] * v2[2] - v1[2] * v2[1];
    r[1] = v1[2] * v2[0] - v1[0] * v2[2];
    r[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }
  
  static inline F
  dot(const F v1[3], const F v2[3]) { return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
  
  // same as cross but assumes points given in inhomogoeneous coordinates
  // as if v[2] = 1
  static inline void
  cross2(const F v1[2], const F v2[2], F r[3])
  {
    r[0] = v1[1] - v2[1];
    r[1] = v2[0] - v1[0];
    r[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }

  // inhomogeneous points with associated tangents to homogeneous line coefficients
  static inline void 
  point_tangent2line(const F p[2], const F tgt[2], F r[3])
  {
    r[0] = -tgt[1]; // normal vector
    r[1] = tgt[0];
    r[2] = p[0]*tgt[1] - p[1]*tgt[0]; // constant term
  }
};

// Not performance critical ---------------------------------------------------
template <typename F>
struct minus_util {
  // Random unit array v of dimension n
  // only on real coordinates, with 0 complex ones
  // we are guaranteeing unifom sampling on the sphere,
  // but simpler rand() on each dimension then normalization also works
  static inline void 
  rand_sphere(C<F> *v/*[chicago14a: 5 minimum, can be 7]*/, unsigned n) {
    F m=0;
    for (unsigned i=0; i < n; ++i) {
      F r = gauss(rnd);
      v[i] = C<F>{r};
      m += r*r;
    }
    // normalization doesn't affect much, just guarantees even sample
    m = std::sqrt(m);
    for (unsigned i=0; i < n; ++i)
      v[i] /= m;
  }

  // random complex
  static void randc(C<F> * __restrict z) { *z = C<F>{gauss(rnd), gauss(rnd)}; *z /= std::abs(*z); }
  static std::random_device rd;
  static std::mt19937 rnd;
  static std::normal_distribution<F> gauss;
  
  // The quaternion order used in the following functions.
  // this is for historical reasons the one used in chicago problem.
  // 
  // Eigen and VXL use {x,y,z,w}, but these functions use {w, x, y z} 
  // Ceres also uses {w, x, y, z}, where w is the scalar part.
  // 
  // Used to cast a q[4] vector to interpret entries and use algorithms not
  // mattering the memory order
  struct quat_shape { F w; F x; F y; F z; };
  
  static inline void normalize_quat(F q[4])
  {
    const F norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    q[0] /= norm; q[1] /= norm; q[2] /= norm; q[3] /= norm;
  }
  
  // arbitrary quaternion to rotation matrix.
  // will normalize the quaternion in-place.
  // based on VXL/VNL
  static inline void quat2rotm(F qq[4], F r[9])
  {
    normalize_quat(qq);
    quat_shape *q = (quat_shape *) qq;
    const F 
      x2 = q->x * q->x,  xy = q->x * q->y,  wx = q->w * q->x,
      y2 = q->y * q->y,  yz = q->y * q->z,  wy = q->w * q->y,
      z2 = q->z * q->z,  zx = q->z * q->x,  wz = q->w * q->z,
      w2 = q->w * q->w;
      
    *r++ = w2 + x2 - y2 - z2;    //  rot(0,0) = r[0]
    *r++ = F(2) * (xy - wz);     //  rot(0,1) = r[1] 
    *r++ = F(2) * (zx + wy);     //  rot(0,2) = r[2] 
    *r++ = F(2) * (xy + wz);     //  rot(1,0) = r[3] 
    *r++ = w2 - x2 + y2 - z2;    //  rot(1,1) = r[4]
    *r++ = F(2) * (yz - wx);     //  rot(1,2) = r[5] 
    *r++ = F(2) * (zx - wy);     //  rot(2,0) = r[6] 
    *r++ = F(2) * (yz + wx);     //  rot(2,1) = r[7] 
    *r   = w2 - x2 - y2 + z2;    //  rot(2,2) = r[8]
  }
  
  // always row-major
  // originally based on Eigen
  static inline void rotm2quat(const F rr[9], F qq[4])
  {
    // use a struct to reinterpret q
    quat_shape *q = (quat_shape *) qq;
    const F (*r)[3] = (const F (*)[3]) rr;
    // coeff_eigen[i] = index in our quaternion shape of corresponding element
    static constexpr unsigned coeff_eigen[4] = {1, 2, 3, 0};

    // This algorithm comes from  "Quaternion Calculus and Fast Animation",
    // Ken Shoemake, 1987 SIGGRAPH course notes
    F t = rr[0] + rr[4] + rr[8]; // trace
    if (t > F(0)) {
      t = std::sqrt(t + F(1.0));
      q->w = F(0.5)*t;
      t = F(0.5)/t;
      q->x = (r[2][1] - r[1][2]) * t;
      q->y = (r[0][2] - r[2][0]) * t;
      q->z = (r[1][0] - r[0][1]) * t;
    } else {
      unsigned i = 0;
      if (r[1][1] > r[0][0]) i = 1;
      if (r[2][2] > r[i][i]) i = 2;
      unsigned j = (i+1)%3;
      unsigned k = (j+1)%3;
      t = std::sqrt(r[i][i]-r[j][j]-r[k][k] + F(1.0));
      qq[coeff_eigen[i]] = F(0.5) * t;
      t = F(0.5)/t;
      q->w = (r[k][j]-r[j][k])*t;
      qq[coeff_eigen[j]] = (r[j][i]+r[i][j])*t;
      qq[coeff_eigen[k]] = (r[k][i]+r[i][k])*t;
    }
  }

  // computes the relative unit quaternion between two unit quaternions
  // based on Eigen quat_product
  // a*conj(b)
  static inline void dquat(const F aa[4], const F bb[4], F d[4])
  {
    const quat_shape *a = (quat_shape *) aa, 
                     *b = (quat_shape *) bb;
    
     *d++ =  a->w * b->w + a->x * b->x + a->y * b->y + a->z * b->z;
     *d++ = -a->w * b->x + a->x * b->w - a->y * b->z + a->z * b->y;
     *d++ = -a->w * b->y + a->y * b->w - a->z * b->x + a->x * b->z;
     *d   = -a->w * b->z + a->z * b->w - a->x * b->y + a->y * b->x;
  }

  // based on Eigen
  static inline void quat_transform(const F q[4], const F v[3], F vrot[3])
  {
    // q*v*q.conj();
    // Note that this algorithm comes from the optimization by hand
    // of the conversion to a Matrix followed by a Matrix/Vector product.
    // It appears to be much faster than the common algorithm found
    // in the literature (30 versus 39 flops). It also requires two
    // Vector3 as temporaries.

    // Vector3 uv = this->vec().cross(v);
    F uv[3];
    minus_3d<F>::cross(q+1, v, uv);
    uv[0] += uv[0]; uv[1] += uv[1]; uv[2] += uv[2];    // uv += uv
    
    // return v + q->w() * uv + q->vec().cross(uv);
    minus_3d<F>::cross(q+1, uv, vrot);
    vrot[0] = v[0] + q[0]*uv[0] + vrot[0];
    vrot[1] = v[1] + q[0]*uv[1] + vrot[1];
    vrot[2] = v[2] + q[0]*uv[2] + vrot[2];
  }
  
  // returns the angle of the smallest rotation around an axis, 
  // such that rotations A and B align. See S. Bianco, G. Ciocca, and D. Marelli,
  // “Evaluating the performance of structure from motion pipelines,” Journal of
  // Imaging, vol. 4, no. 8, 2018
  //
  // Input: unit quaternions qa and qb
  //
  //rotation_error(const F Ra[ncoords3d][ncoords3d], const F Rb[ncoords3d][ncoords3d])
  //{
  //    dR = norm(skew2v(Rots{n}*R_tilde'));
  //}
  // Based on Eigen
  static inline F rotation_error(const F p[4], const F q[4])
  {
    // normalize_quat(p); normalize_quat(q);
    F d[4];
    dquat(p, q , d);

    const F vnorm = sqrt(d[1]*d[1] + d[2]*d[2] + d[3]*d[3]);
    return F(2) * std::atan2(vnorm, std::fabs(d[0]));
  }
};

} // namespace minus

#endif // internal_util_h_
