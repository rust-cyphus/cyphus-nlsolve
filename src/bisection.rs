enum NLSolveCode {
    Success,
    Continue,
    InvalidArgument,
    BadFunction,
}

pub struct Bisection<F: Fn(f64) -> f64> {
    f: F,
    lower: f64,
    upper: f64,
    f_lower: f64,
    f_upper: f64,
}

impl<F: Fn(f64) -> f64> Bisection<F> {
    fn new(f: F, lower: f64, upper: f64) -> Result<Bisection<F>, NLSolveCode> {
        let root = 0.5 * (lower + upper);
        let f_upper = f(upper);
        let f_lower = f(lower);

        if f_upper.is_nan() || f_lower.is_nan() {
            Err(NLSolveCode::BadFunction)
        } else if f_lower * f_upper > 0.0 {
            Err(NLSolveCode::InvalidArgument)
        } else {
            Ok(Bisection {
                f,
                lower,
                upper,
                f_lower,
                f_upper,
            })
        }
    }

    fn iterate(&mut self) -> NLSolveCode {
        let right = self.upper;
        let left = self.lower;

        if self.f_lower == 0.0 {
            self.upper = left;
            return NLSolveCode::Success;
        }
        if self.f_upper == 0.0 {
            self.lower = right;
            return NLSolveCode::Success;
        }

        let bisect = (left + right) / 2.0;
        let fbisect = (self.f)(bisect);

        if fbisect.is_nan() {
            return NLSolveCode::BadFunction;
        }

        if fbisect == 0.0 {
            self.lower = bisect;
            self.upper = bisect;
            return NLSolveCode::Success;
        }

        // Discard the half of the interval which doesn't contian the root
        if self.f_lower * fbisect < 0.0 {
            self.upper = bisect;
            self.f_upper = fbisect;
        } else {
            self.lower = bisect;
            self.f_lower = fbisect;
        }
        return NLSolveCode::Success;
    }

    pub fn solve(f: F, a: f64, b: f64, tol: f64) -> f64 {
        let mut solver = Bisection::new(&f, a, b);

        match solver {
            Ok(ref mut sol) => loop {
                sol.iterate();
                if sol.f_lower.abs() <= tol || sol.f_upper.abs() <= tol {
                    // linearly interpolate
                    return (sol.upper * sol.f_lower - sol.lower * sol.f_upper)
                        / (sol.f_lower - sol.f_upper);
                }
            },
            _ => {
                // Need to issue warning
                return (a + b) / 2.0;
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_trig() {
        let tol = 1e-5;
        let mut root = Bisection::solve(|x| x.sin(), 3.0, 4.0, tol);
        assert!((root - std::f64::consts::PI).abs() < tol);

        root = Bisection::solve(|x| x.sin(), -4.0, -3.0, tol);
        assert!((root + std::f64::consts::PI).abs() < tol);

        root = Bisection::solve(|x| x.sin(), -1.0 / 3.0, 1.0, tol);
        assert!(root.abs() < tol);

        root = Bisection::solve(|x| x.cos(), 0.0, 3.0, tol);
        assert!((root - std::f64::consts::PI / 2.0).abs() < tol);

        root = Bisection::solve(|x| x.cos(), -3.0, 0.0, tol);
        assert!((root + std::f64::consts::PI / 2.0).abs() < tol);
    }
}
