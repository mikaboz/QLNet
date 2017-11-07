using System;
using System.Collections.Generic;
using System.Linq;

namespace QLNet
{
   public class CIR2 : TwoFactorAffineModel
   {

      public class FellerConstraint : Constraint
      {
         private class Impl : IConstraint
         {
            public bool test(Vector param)
            {
               double kappa1 = param[0];
               double theta1 = param[1];
               double sigma1 = param[2];

               double kappa2 = param[3];
               double theta2 = param[4];
               double sigma2 = param[5];


               return (sigma1 >= 0.0 && sigma1 * sigma1 < 2.0 * kappa1 * theta1) || (sigma2 >= 0.0 && sigma2 * sigma2 < 2.0 * kappa2 * theta2);
            }

            public Vector upperBound(Vector parameters)
            {
               return new Vector(parameters.size(), Double.MaxValue);
            }

            public Vector lowerBound(Vector parameters)
            {
               return new Vector(parameters.size(), Double.MinValue);
            }
         }

         public FellerConstraint()
            : base(new FellerConstraint.Impl())
         { }

      }

      #region Accessors
      public double X0 { get { return First.r0_; } }
      public double Y0 { get { return Second.r0_; } }
      public new CoxIngersollRoss First { get { return (CoxIngersollRoss)base.First; } }
      public new CoxIngersollRoss Second { get { return (CoxIngersollRoss)base.Second; } }
      public double Kappa1 { get { return this.First.Kappa; } }
      public double Theta1 {get { return this.First.Theta; } }
      public double Sigma1 { get { return this.First.Sigma; } }
      public double Kappa2 { get { return this.Second.Kappa; } }
      public double Theta2 { get { return this.Second.Theta; } }
      public double Sigma2 { get { return this.Second.Sigma; } }
      #endregion

      #region Constructors
      public CIR2(double x0, double kappa1, double theta1, double sigma1, double y0, double kappa2, double theta2, double sigma2) :
         this(new CoxIngersollRoss(x0,kappa1,theta1,sigma1), new CoxIngersollRoss(y0, kappa2, theta2, sigma2))
      {}
      public CIR2(CoxIngersollRoss first, CoxIngersollRoss second) :
         base(first, second,null)
      {}
      #endregion

      #region IAffineModel implémentation
      public override double A(double t, double T)
      {
         return First.A(t, T) * Second.A(t, T);
      }
      #endregion

      /// <summary>
      /// TODO Cf. Brigo&Mercurio p.177
      /// </summary>
      public override double DiscountBondOption(Option.Type type, double strike, double maturity, double bondMaturity)
      {
         throw new NotImplementedException();
      }
   }
}
