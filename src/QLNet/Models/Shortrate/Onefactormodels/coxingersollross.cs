/*
 Copyright (C) 2010 Philippe Real (ph_real@hotmail.com)
  
 This file is part of QLNet Project https://github.com/amaggiulli/qlnet

 QLNet is free software: you can redistribute it and/or modify it
 under the terms of the QLNet license.  You should have received a
 copy of the license along with this program; if not, license is  
 available online at <http://qlnet.sourceforge.net/License.html>.
  
 QLNet is a based on QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/
 The QuantLib license is available online at http://quantlib.org/license.shtml.
 
 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/
using System;

namespace QLNet
{
   /// <summary>
   /// Cox-Ingersoll-Ross model class.
   /// <remarks>
   /// This class implements the Cox-Ingersoll-Ross model defined by
   /// dr_t = k(\theta - r_t)dt + \sqrt{r_t}\sigma dW_t .
   /// This class was not tested enough to guarantee its functionality.
   /// </remarks>
   /// </summary>
   public class CoxIngersollRoss : OneFactorAffineModel
   {
      #region Feller constraint
      public class FellerConstraint : Constraint
      {
         private class Impl : IConstraint
         {
            public bool test(Vector param)
            {
               double theta = param[0];
               double kappa = param[1];
               double sigma = param[2];

               return (sigma >= 0.0 && sigma * sigma < 2.0 * kappa * theta);
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
      #endregion

      #region Constructors
      public CoxIngersollRoss(double r0, double kappa = 0.1, double theta = 0.1, double sigma = 0.1) :
         base(3)
      {
         Utils.QL_REQUIRE(r0 >= 0, () => "r0 must be positive to initially satisfy feller constraint");
         constraint_ = new CompositeConstraint(base.constraint_, new FellerConstraint());
         r0_ = r0;
         arguments_[0] = new ConstantParameter(kappa, new PositiveConstraint());
         arguments_[1] = new ConstantParameter(theta, new PositiveConstraint());
         arguments_[2] = new ConstantParameter(sigma, new PositiveConstraint());
      }
      public CoxIngersollRoss(Handle<YieldTermStructure> termStructure, double kappa = 0.1, double theta = 0.1, double sigma = 0.1) :
         this(termStructure.link.forwardRate(0,0,Compounding.Continuous,Frequency.NoFrequency).rate(),kappa,theta,sigma)
      {}
      #endregion

      #region Accessors
      public double r0_;
      public double Kappa
      {
         get { return arguments_[0].value(0.0); }
      }
      public double Theta
      {
         get { return arguments_[1].value(0.0); }
      }
      public double Sigma
      {
         get { return arguments_[2].value(0.0); }
      }
      #endregion

      #region IAffineModel implémentation
      public override double DiscountBondOption(Option.Type type, double strike, double maturity, double bondMaturity)
      {
         Utils.QL_REQUIRE(strike > 0.0, () => "strike must be positive");
         double discountT = this.DiscountBond(0.0, maturity, r0_);
         double discountS = this.DiscountBond(0.0, bondMaturity, r0_);

         if (maturity < Const.QL_EPSILON)
         {
            switch (type)
            {
               case Option.Type.Call:
                  return Math.Max(discountS - strike, 0.0);
               case Option.Type.Put:
                  return Math.Max(strike - discountS, 0.0);
               default:
                  Utils.QL_FAIL("unsupported option type");
                  break;
            }
         }
         double sigma2 = Sigma * Sigma;
         double h = Math.Sqrt(Kappa * Kappa + 2.0 * sigma2);
         double b = B(maturity, bondMaturity);

         double rho = 2.0 * h / (sigma2 * (Math.Exp(h * maturity) - 1.0));
         double psi = (Kappa + h) / sigma2;

         double df = 4.0 * Kappa * Theta / sigma2;
         double ncps = 2.0 * rho * rho * r0_ * Math.Exp(h * maturity) / (rho + psi + b);
         double ncpt = 2.0 * rho * rho * r0_ * Math.Exp(h * maturity) / (rho + psi);

         NonCentralChiSquareDistribution chis = new NonCentralChiSquareDistribution(df, ncps);
         NonCentralChiSquareDistribution chit = new NonCentralChiSquareDistribution(df, ncpt);

         double z = Math.Log(A(maturity, bondMaturity) / strike) / b;
         double call = discountS * chis.value(2.0 * z * (rho + psi + b)) -
                       strike * discountT * chit.value(2.0 * z * (rho + psi));

         if (type == Option.Type.Call)
            return call;
         else
            return call - discountS + strike * discountT;
      }
      public override double A(double t, double T)
      {
         double sigma2 = Sigma * Sigma;
         double h = Math.Sqrt(Kappa * Kappa + 2.0 * sigma2);
         double numerator = 2.0 * h * Math.Exp(0.5 * (Kappa + h) * (T - t));
         double denominator = 2.0 * h + (Kappa + h) * (Math.Exp((T - t) * h) - 1.0);
         double value = Math.Log(numerator / denominator) *
                        2.0 * Kappa * Theta / sigma2;
         return Math.Exp(value);
      }
      public override double B(double t, double T)
      {
         double h = Math.Sqrt(Kappa * Kappa + 2.0 * Sigma * Sigma);
         double temp = Math.Exp((T - t) * h) - 1.0;
         double numerator = 2.0 * temp;
         double denominator = 2.0 * h + (Kappa + h) * temp;
         double value = numerator / denominator;
         return value;
      }
      #endregion

      #region Dynamics
      public override ShortRateModel.Dynamics dynamics()
      {
         return new OneFactorModel.Dynamics(new SquareRootProcess(r0_, Kappa, Theta, Sigma));
      }
      
      #endregion

      #region Tree
      public override Lattice tree(TimeGrid grid)
      {
         OneFactorModel.Dynamics dyn = (OneFactorModel.Dynamics)dynamics();
         TrinomialTree trinomial = new TrinomialTree(dyn.Process, grid, true);
         return new ShortRateTree(trinomial, dyn, grid);
      }
      #endregion
   }
}
