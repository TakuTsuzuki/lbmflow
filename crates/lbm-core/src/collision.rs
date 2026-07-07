//! Collision composition boundary for the BGK/TRT + Guo path.
//!
//! This is the first mechanical seam: the operator is written against an
//! arithmetic sink, while the central-moment branch remains on its existing
//! direct implementation.

use crate::kernels::RawSlice;
use crate::lattice::{Lattice, Q_MAX};
use crate::params::KParams;
use crate::real::Real;

pub(crate) trait Arith {
    type V: Copy;

    fn lit(&mut self, x: f64) -> Self::V;
    fn add(&mut self, a: Self::V, b: Self::V) -> Self::V;
    fn sub(&mut self, a: Self::V, b: Self::V) -> Self::V;
    fn mul(&mut self, a: Self::V, b: Self::V) -> Self::V;

    fn pop(&mut self, q: usize) -> Self::V;
    fn set_pop(&mut self, q: usize, v: Self::V);
    fn rho(&mut self) -> Self::V;
    fn u(&mut self, d: usize) -> Self::V;
    fn force(&mut self, d: usize) -> Self::V;
    fn omega_p(&mut self) -> Self::V;
    fn omega_m(&mut self) -> Self::V;
    fn cp(&mut self) -> Self::V;
    fn cm(&mut self) -> Self::V;
    fn force_on(&self) -> bool;
    fn w(&self, q: usize) -> f64;
    fn c(&self, q: usize, d: usize) -> i8;
}

pub(crate) trait Equilibrium {
    fn feq<A: Arith, L: Lattice>(a: &mut A, q: usize, usq: A::V) -> A::V;
}

pub(crate) trait Forcing {
    fn source<A: Arith, L: Lattice>(a: &mut A, q: usize, uf: A::V, cu: A::V) -> A::V;
}

pub(crate) trait Collision {
    type Eq: Equilibrium;
    type Force: Forcing;

    fn relax<A: Arith, L: Lattice>(a: &mut A);
}

pub(crate) struct SecondOrderEq;
pub(crate) struct GuoForcing;
pub(crate) struct TrtGuo;

#[inline(always)]
fn dot_u<A: Arith, L: Lattice>(a: &mut A, q: usize) -> A::V {
    let c0 = a.lit(a.c(q, 0) as f64);
    let u0 = a.u(0);
    let mut acc = a.mul(c0, u0);
    for d in 1..L::D {
        let cd = a.lit(a.c(q, d) as f64);
        let ud = a.u(d);
        let term = a.mul(cd, ud);
        acc = a.add(acc, term);
    }
    acc
}

#[inline(always)]
fn dot_f<A: Arith, L: Lattice>(a: &mut A, q: usize) -> A::V {
    let c0 = a.lit(a.c(q, 0) as f64);
    let f0 = a.force(0);
    let mut acc = a.mul(c0, f0);
    for d in 1..L::D {
        let cd = a.lit(a.c(q, d) as f64);
        let fd = a.force(d);
        let term = a.mul(cd, fd);
        acc = a.add(acc, term);
    }
    acc
}

impl Equilibrium for SecondOrderEq {
    #[inline(always)]
    fn feq<A: Arith, L: Lattice>(a: &mut A, q: usize, usq: A::V) -> A::V {
        let cu = dot_u::<A, L>(a, q);
        let rho = a.rho();
        let one = a.lit(1.0);
        let drho = a.sub(rho, one);
        let three = a.lit(3.0);
        let three_cu = a.mul(three, cu);
        let f45 = a.lit(4.5);
        let f45_cu = a.mul(f45, cu);
        let f45_cu2 = a.mul(f45_cu, cu);
        let f15 = a.lit(1.5);
        let f15_usq = a.mul(f15, usq);
        let sum = a.add(three_cu, f45_cu2);
        let inner = a.sub(sum, f15_usq);
        let rho_inner = a.mul(rho, inner);
        let body = a.add(drho, rho_inner);
        let w = a.lit(a.w(q));
        a.mul(w, body)
    }
}

impl Forcing for GuoForcing {
    #[inline(always)]
    fn source<A: Arith, L: Lattice>(a: &mut A, q: usize, uf: A::V, cu: A::V) -> A::V {
        if !a.force_on() {
            return a.lit(0.0);
        }
        let cf = dot_f::<A, L>(a, q);
        let cf_minus_uf = a.sub(cf, uf);
        let three = a.lit(3.0);
        let first = a.mul(three, cf_minus_uf);
        let nine = a.lit(9.0);
        let nine_cu = a.mul(nine, cu);
        let second = a.mul(nine_cu, cf);
        let body = a.add(first, second);
        let w = a.lit(a.w(q));
        a.mul(w, body)
    }
}

impl Collision for TrtGuo {
    type Eq = SecondOrderEq;
    type Force = GuoForcing;

    #[inline(always)]
    fn relax<A: Arith, L: Lattice>(a: &mut A) {
        let half = a.lit(0.5);
        let u0 = a.u(0);
        let mut usq = a.mul(u0, u0);
        for d in 1..L::D {
            let ud = a.u(d);
            let ud2 = a.mul(ud, ud);
            usq = a.add(usq, ud2);
        }
        let f0 = a.force(0);
        let mut uf = a.mul(u0, f0);
        for d in 1..L::D {
            let ud = a.u(d);
            let fd = a.force(d);
            let udfd = a.mul(ud, fd);
            uf = a.add(uf, udfd);
        }

        let mut feq = [a.lit(0.0); Q_MAX];
        let mut src = [a.lit(0.0); Q_MAX];
        for q in 0..L::Q {
            let cu = dot_u::<A, L>(a, q);
            feq[q] = Self::Eq::feq::<A, L>(a, q, usq);
            src[q] = Self::Force::source::<A, L>(a, q, uf, cu);
        }

        let rest = L::REST;
        let f0 = a.pop(rest);
        let op = a.omega_p();
        let cp = a.cp();
        let f0_minus_e = a.sub(f0, feq[rest]);
        let op_delta = a.mul(op, f0_minus_e);
        let relaxed = a.sub(f0, op_delta);
        let cp_src = a.mul(cp, src[rest]);
        let rest_out = a.add(relaxed, cp_src);
        a.set_pop(rest, rest_out);

        for &(qa, qb) in L::PAIRS {
            let fa = a.pop(qa);
            let fb = a.pop(qb);
            let fa_plus_fb = a.add(fa, fb);
            let fp = a.mul(half, fa_plus_fb);
            let fa_minus_fb = a.sub(fa, fb);
            let fm = a.mul(half, fa_minus_fb);
            let e_plus = a.add(feq[qa], feq[qb]);
            let ep = a.mul(half, e_plus);
            let e_minus = a.sub(feq[qa], feq[qb]);
            let em = a.mul(half, e_minus);
            let s_plus = a.add(src[qa], src[qb]);
            let sp = a.mul(half, s_plus);
            let s_minus = a.sub(src[qa], src[qb]);
            let sm = a.mul(half, s_minus);
            let fp_minus_ep = a.sub(fp, ep);
            let op = a.omega_p();
            let rp = a.mul(op, fp_minus_ep);
            let fm_minus_em = a.sub(fm, em);
            let om = a.omega_m();
            let rm = a.mul(om, fm_minus_em);
            let cp = a.cp();
            let cp_sp = a.mul(cp, sp);
            let cm = a.cm();
            let cm_sm = a.mul(cm, sm);
            let fa_rp = a.sub(fa, rp);
            let fa_rp_rm = a.sub(fa_rp, rm);
            let src_a = a.add(cp_sp, cm_sm);
            let va = a.add(fa_rp_rm, src_a);
            let fb_rp = a.sub(fb, rp);
            let fb_rp_rm = a.add(fb_rp, rm);
            let src_b = a.sub(cp_sp, cm_sm);
            let vb = a.add(fb_rp_rm, src_b);
            a.set_pop(qa, va);
            a.set_pop(qb, vb);
        }
    }
}

pub(crate) struct ScalarArith<'a, T: Real> {
    f: RawSlice<T>,
    np: usize,
    i: usize,
    rho: T,
    u: [T; 3],
    force: [T; 3],
    omega_p: T,
    omega_m: T,
    cp: T,
    cm: T,
    force_on: bool,
    kp: &'a KParams<T>,
}

impl<'a, T: Real> ScalarArith<'a, T> {
    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    pub(crate) fn new(
        f: RawSlice<T>,
        np: usize,
        i: usize,
        rho: T,
        u: [T; 3],
        force: [T; 3],
        omega_p: T,
        cp: T,
        force_on: bool,
        kp: &'a KParams<T>,
    ) -> Self {
        Self {
            f,
            np,
            i,
            rho,
            u,
            force,
            omega_p,
            omega_m: kp.omega_m,
            cp,
            cm: kp.cm,
            force_on,
            kp,
        }
    }
}

impl<T: Real> Arith for ScalarArith<'_, T> {
    type V = T;

    #[inline(always)]
    fn lit(&mut self, x: f64) -> Self::V {
        T::r(x)
    }

    #[inline(always)]
    fn add(&mut self, a: Self::V, b: Self::V) -> Self::V {
        a + b
    }

    #[inline(always)]
    fn sub(&mut self, a: Self::V, b: Self::V) -> Self::V {
        a - b
    }

    #[inline(always)]
    fn mul(&mut self, a: Self::V, b: Self::V) -> Self::V {
        a * b
    }

    #[inline(always)]
    fn pop(&mut self, q: usize) -> Self::V {
        // SAFETY: caller provides a row-disjoint cell under the RawSlice contract.
        unsafe { self.f.get(q * self.np + self.i) }
    }

    #[inline(always)]
    fn set_pop(&mut self, q: usize, v: Self::V) {
        // SAFETY: caller provides a row-disjoint cell under the RawSlice contract.
        unsafe { self.f.set(q * self.np + self.i, v) };
    }

    #[inline(always)]
    fn rho(&mut self) -> Self::V {
        self.rho
    }

    #[inline(always)]
    fn u(&mut self, d: usize) -> Self::V {
        self.u[d]
    }

    #[inline(always)]
    fn force(&mut self, d: usize) -> Self::V {
        self.force[d]
    }

    #[inline(always)]
    fn omega_p(&mut self) -> Self::V {
        self.omega_p
    }

    #[inline(always)]
    fn omega_m(&mut self) -> Self::V {
        self.omega_m
    }

    #[inline(always)]
    fn cp(&mut self) -> Self::V {
        self.cp
    }

    #[inline(always)]
    fn cm(&mut self) -> Self::V {
        self.cm
    }

    #[inline(always)]
    fn force_on(&self) -> bool {
        self.force_on
    }

    #[inline(always)]
    fn w(&self, q: usize) -> f64 {
        self.kp.wr[q].as_f64()
    }

    #[inline(always)]
    fn c(&self, q: usize, d: usize) -> i8 {
        self.kp.cr[q][d].as_f64() as i8
    }
}
