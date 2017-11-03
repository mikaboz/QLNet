using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QLNet
{
   /*
   public class GeneralizedGaussianModel : MultiFactorAffineModel, ITermStructureConsistentModel
   {
      #region ITermStructureConsistentModel
      protected Handle<YieldTermStructure> termStructure_;
      public Handle<YieldTermStructure> TermStructure { get { return termStructure_; } }
      #endregion

      #region Fitting
      FittingParameter phi_;
      protected override void generateArguments()
      {
         phi_ = new FittingParameter(this);
      }
      private double ValueForFitting(double t)
      {
         double doubleSum = 0;
         for (int i = 0; i < nFactors; i++)
            for (int j = 0; j < nFactors; j++)
            {
               double expi = 1 - Math.Exp(-Kappa(i) * t);
               double expj = 1 - Math.Exp(-Kappa(j) * t);
               doubleSum += this.Rho(i, j) * Sigma(i) * Sigma(j) / (Kappa(i) * Kappa(j)) * expi * expj;
            }
         return 0.5 * doubleSum;
      }
      public class FittingParameter : TermStructureFittingParameter
      {
         private new class Impl : TermStructureFittingParameter.Impl
         {
            public Impl(GeneralizedGaussianModel model) :
               base(model)
            { }

            public override double value(Vector v, double t)
            {
               double forward = model_.TermStructureInitialForwardRate(t);

               GeneralizedGaussianModel g2 = model_ as GeneralizedGaussianModel;
               return forward + g2.ValueForFitting(t);
            }
         }
         public FittingParameter(GeneralizedGaussianModel model) :
            base(new Impl(model))
         { }
      }
      #endregion

      #region Accessors
      Gaussian Factor(int i) { return factors_[i] as Gaussian; }
      double Kappa(int i) { return Factor(i).Arguments[0].value(0.0); }
      double Sigma(int i) { return Factor(i).Arguments[1].value(0.0); }
      #endregion

      #region Constructors
      public static class Generators
      {
         public static IEnumerable<OneFactorAffineModel> GenerateFactors(double[] kappa, double[] sigma)
         {
            Utils.QL_REQUIRE(kappa.Length == sigma.Length, () => "kappa et sigma doivent être de meme tailles.");
            List<OneFactorAffineModel> list = new List<OneFactorAffineModel>();
            for (int i = 0; i < kappa.Length; i++)
               list.Add(new Gaussian(kappa[i], sigma[i]));
            return list.AsEnumerable();
         }
         public static IEnumerable<OneFactorAffineModel> GenerateFactors(int nFactors)
         {
            double[] kappa = new double[nFactors];
            double[] sigma = new double[nFactors];
            for (int i = 0; i < nFactors; i++)
            {
               kappa[i] = 0.1;
               sigma[i] = 0.2;
            }
            return GenerateFactors(kappa, sigma);

         }
         public static double[,] GenerateWeakCorrelation(int nFactors)
         {
            double[,] cor = new double[nFactors, nFactors];
            for (int i = 0; i < nFactors; i++)
               cor[i, i] = 1.0;
            return cor;
         }

      }
      public GeneralizedGaussianModel(Handle<YieldTermStructure> termStructure, double[] kappa, double[] sigma, double[,] rho) :
         base(Generators.GenerateFactors(kappa,sigma),rho)
      {
         termStructure_ = termStructure;
         termStructure.registerWith(update);
         generateArguments();
      }
      /// <summary>
      /// Gn++ model initialized with default values kappa_i = 0.1 , sigma_i = 0.2 , rho_i_j = 0
      /// </summary>
      public GeneralizedGaussianModel(Handle<YieldTermStructure> termStructure, int nFactors) :
         base(Generators.GenerateFactors(nFactors), Generators.GenerateWeakCorrelation(nFactors))
      {
         termStructure_ = termStructure;
         termStructure.registerWith(update);
         generateArguments();
      }
      #endregion

      #region Dynamics
      public override ShortRateModel.Dynamics dynamics()
      {
         return new GeneralizedGaussianModel.Dynamics(this);
      }
      public new class Dynamics : MultiFactorModel.Dynamics
      {
         Parameter fitting_;
         public Dynamics(GeneralizedGaussianModel model) :
            base(model)
         {
            fitting_ = model.phi_;
         }
         public override double shortRate(double t, Vector variables)
         {
            return fitting_.value(t) + base.shortRate(t, variables);
         }
      }
      #endregion

      #region IAffineModel override
      private double V(double t, double T)
      {
         double doubleSum = 0;
         double dt = T - t;
         for (int i = 0; i < nFactors; i++)
            for (int j = 0; j < nFactors; j++)
            {
               double expi = (1 - Math.Exp(-Kappa(i) * dt)) / Kappa(i);
               double expj = (1 - Math.Exp(-Kappa(j) * dt)) / Kappa(j);
               double expij = (1 - Math.Exp(-(Kappa(i) + Kappa(j)) * dt)) / (Kappa(i) + Kappa(j));
               doubleSum += this.Rho(i, j) * Sigma(i) * Sigma(j) / (Kappa(i) * Kappa(j)) * (dt - expi - expj + expij);
            }
         return doubleSum;
      }
      public override double A(double t, double T)
      {
         double PM0t = termStructure_.link.discount(t);
         double PM0T = termStructure_.link.discount(T);
         return PM0T / PM0t * Math.Exp(0.5 * (V(t, T) - V(0, T) + V(0, t)));
      }


      public virtual double discountBondOption(Option.Type type, double strike, double maturity, double bondMaturity)
      {
         double v = sigmaP(maturity, bondMaturity);
         double f = termStructure_.link.discount(bondMaturity);
         double k = termStructure_.link.discount(maturity) * strike;

         return Utils.blackFormula(type, k, f, v);
      }
      #endregion

      double sigmaP(double t, double s)
      {
         double temp = 1.0 - Math.Exp(-(Kappa1 + Kappa2) * t);
         double temp1 = 1.0 - Math.Exp(-Kappa1 * (s - t));
         double temp2 = 1.0 - Math.Exp(-Kappa2 * (s - t));
         double a3 = Kappa1 * Kappa1 * Kappa1;
         double b3 = Kappa2 * Kappa2 * Kappa2;
         double sigma2 = Sigma1 * Sigma1;
         double eta2 = Sigma2 * Sigma2;
         double value =
             0.5 * sigma2 * temp1 * temp1 * (1.0 - Math.Exp(-2.0 * Kappa2 * t)) / a3 +
             0.5 * eta2 * temp2 * temp2 * (1.0 - Math.Exp(-2.0 * Kappa2 * t)) / b3 +
             2.0 * Rho * Sigma1 * Sigma2 / (Kappa1 * Kappa2 * (Kappa1 + Kappa2)) *
             temp1 * temp2 * temp;
         return Math.Sqrt(value);
      }

      public double swaption(Swaption.Arguments arguments, double fixedRate, double range, int intervals)
      {
         Date settlement = termStructure_.link.referenceDate();
         DayCounter dayCounter = termStructure_.link.dayCounter();
         double start = dayCounter.yearFraction(settlement,
                                              arguments.floatingResetDates[0]);
         double w = (arguments.type == VanillaSwap.Type.Payer ? 1 : -1);

         List<double> fixedPayTimes = new InitializedList<double>(arguments.fixedPayDates.Count);
         for (int i = 0; i < fixedPayTimes.Count; ++i)
            fixedPayTimes[i] =
                dayCounter.yearFraction(settlement,
                                        arguments.fixedPayDates[i]);

         SwaptionPricingFunction function = new SwaptionPricingFunction(a(),
                                                 sigma(), b(), eta(), rho(),
                                                 w, start,
                                                 fixedPayTimes,
                                                 fixedRate, this);

         double upper = function.mux() + range * function.sigmax();
         double lower = function.mux() - range * function.sigmax();
         SegmentIntegral integrator = new SegmentIntegral(intervals);
         return arguments.nominal * w * termStructure_.link.discount(start) *
             integrator.value(function.value, lower, upper);
      }

      public class SwaptionPricingFunction
      {

         #region private fields 
         double a_, sigma_, b_, eta_, rho_, w_;
         double T_;
         List<double> t_;
         double rate_;
         int size_;
         Vector A_, Ba_, Bb_;
         double mux_, muy_, sigmax_, sigmay_, rhoxy_;
         #endregion

         public SwaptionPricingFunction(double a, double sigma,
                                 double b, double eta, double rho,
                                 double w, double start,
                                 List<double> payTimes,
                                 double fixedRate, G2 model)
         {
            a_ = a;
            sigma_ = sigma;
            b_ = b;
            eta_ = eta;
            rho_ = rho;
            w_ = w;
            T_ = start;
            t_ = payTimes;
            rate_ = fixedRate;
            size_ = t_.Count;

            A_ = new Vector(size_);
            Ba_ = new Vector(size_);
            Bb_ = new Vector(size_);


            sigmax_ = sigma_ * Math.Sqrt(0.5 * (1.0 - Math.Exp(-2.0 * a_ * T_)) / a_);
            sigmay_ = eta_ * Math.Sqrt(0.5 * (1.0 - Math.Exp(-2.0 * b_ * T_)) / b_);
            rhoxy_ = rho_ * eta_ * sigma_ * (1.0 - Math.Exp(-(a_ + b_) * T_)) /
                ((a_ + b_) * sigmax_ * sigmay_);

            double temp = sigma_ * sigma_ / (a_ * a_);
            mux_ = -((temp + rho_ * sigma_ * eta_ / (a_ * b_)) * (1.0 - Math.Exp(-a * T_)) -
                     0.5 * temp * (1.0 - Math.Exp(-2.0 * a_ * T_)) -
                     rho_ * sigma_ * eta_ / (b_ * (a_ + b_)) *
                     (1.0 - Math.Exp(-(b_ + a_) * T_)));

            temp = eta_ * eta_ / (b_ * b_);
            muy_ = -((temp + rho_ * sigma_ * eta_ / (a_ * b_)) * (1.0 - Math.Exp(-b * T_)) -
                     0.5 * temp * (1.0 - Math.Exp(-2.0 * b_ * T_)) -
                     rho_ * sigma_ * eta_ / (a_ * (a_ + b_)) *
                     (1.0 - Math.Exp(-(b_ + a_) * T_)));

            for (int i = 0; i < size_; i++)
            {
               A_[i] = model.A(T_, t_[i]);
               Ba_[i] = model.B(a_, t_[i] - T_);
               Bb_[i] = model.B(b_, t_[i] - T_);
            }
         }

         internal double mux() { return mux_; }

         internal double sigmax() { return sigmax_; }

         public double value(double x)
         {
            CumulativeNormalDistribution phi = new CumulativeNormalDistribution();
            double temp = (x - mux_) / sigmax_;
            double txy = Math.Sqrt(1.0 - rhoxy_ * rhoxy_);

            Vector lambda = new Vector(size_);
            int i;
            for (i = 0; i < size_; i++)
            {
               double tau = (i == 0 ? t_[0] - T_ : t_[i] - t_[i - 1]);
               double c = (i == size_ - 1 ? (1.0 + rate_ * tau) : rate_ * tau);
               lambda[i] = c * A_[i] * Math.Exp(-Ba_[i] * x);
            }

            SolvingFunction function = new SolvingFunction(lambda, Bb_);
            Brent s1d = new Brent();
            s1d.setMaxEvaluations(1000);
            double yb = s1d.solve(function, 1e-6, 0.00, -100.0, 100.0);

            double h1 = (yb - muy_) / (sigmay_ * txy) -
                rhoxy_ * (x - mux_) / (sigmax_ * txy);
            double value = phi.value(-w_ * h1);


            for (i = 0; i < size_; i++)
            {
               double h2 = h1 +
                   Bb_[i] * sigmay_ * Math.Sqrt(1.0 - rhoxy_ * rhoxy_);
               double kappa = -Bb_[i] *
                   (muy_ - 0.5 * txy * txy * sigmay_ * sigmay_ * Bb_[i] +
                    rhoxy_ * sigmay_ * (x - mux_) / sigmax_);
               value -= lambda[i] * Math.Exp(kappa) * phi.value(-w_ * h2);
            }

            return Math.Exp(-0.5 * temp * temp) * value /
                (sigmax_ * Math.Sqrt(2.0 * QLNet.Const.M_PI));
         }

         public class SolvingFunction : ISolver1d
         {

            Vector lambda_;
            Vector Bb_;

            public SolvingFunction(Vector lambda, Vector Bb)
            {
               lambda_ = lambda;
               Bb_ = Bb;
            }

            public override double value(double y)
            {
               double value = 1.0;
               for (int i = 0; i < lambda_.size(); i++)
               {
                  value -= lambda_[i] * Math.Exp(-Bb_[i] * y);
               }
               return value;
            }

         }

      }

   }
   */
}
