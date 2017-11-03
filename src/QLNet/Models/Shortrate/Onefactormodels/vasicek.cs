/*
 Copyright (C) 2010 Philippe Real (ph_real@hotmail.com)
 Copyright (C) 2008-2016 Andrea Maggiulli (a.maggiulli@gmail.com)
  
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
   /// Vasicek model class
   /// </summary>
   public class Vasicek : OneFactorAffineModel
   {
      protected double r0_;
      public Vasicek(double r0, double kappa, double theta = 0.05, double sigma = 0.01)
         : base(3)
      {
         r0_ = r0;
         arguments_[0] = new ConstantParameter(kappa, new PositiveConstraint());
         arguments_[1] = new ConstantParameter(theta, new NoConstraint());
         arguments_[2] = new ConstantParameter(sigma, new PositiveConstraint());
      }
      public Vasicek(Handle<YieldTermStructure> termStructure, double kappa, double theta = 0.05, double sigma = 0.01)
         : this(termStructure.link.forwardRate(0, 0, Compounding.Continuous, Frequency.NoFrequency).rate(), kappa, theta, sigma)
      { }

      public virtual double Kappa
      {
         get { return arguments_[0].value(0.0); }
      }
      public virtual double Theta
      {
         get { return arguments_[1].value(0.0); }
      }
      public virtual double Sigma
      {
         get { return arguments_[2].value(0.0); }
      }

      protected double V(double t, double T)
      {
         double exp = Math.Exp(-Kappa * (T - t));
         double temp1 = (1 - exp) / Kappa;
         double temp2 = (1-exp*exp)/ Kappa;
         double c = Sigma / Kappa;
         return c * c * (T - t - 2 * temp1 + 0.5 * temp2);
      }

      protected double E1(double t, double T)
      {
         double exp = Math.Exp(-Kappa * (T - t));
         return (1 - exp) / Kappa;
      }
      protected double E2(double t, double T)
      {
         return Theta * (T - t - E1(t, T));
      }

      public override double A(double t, double T)
      {
         return Math.Exp(0.5*V(t, T) - E2(t, T));
         /*
         double _a = Kappa;
         if (_a < Math.Sqrt(Const.QL_EPSILON))
         {
            return 0.0;
         }
         else
         {
            double sigma2 = Sigma * Sigma;
            double bt = B(t, T);
            return Math.Exp((Theta - 0.5 * sigma2 / (_a * _a)) * (bt - (T - t))
                            - 0.25 * sigma2 * bt * bt / _a);
         }
         */
      }
      public override double B(double t, double T)
      {
         return E1(t, T);
         /*
         double _a = Kappa;
         if (_a < Math.Sqrt(Const.QL_EPSILON))
            return (T - t);
         else
            return (1.0 - Math.Exp(-_a * (T - t))) / _a;
         */
      }

      public override double DiscountBondOption(Option.Type type, double strike, double maturity, double bondMaturity)
      {
         double v;
         if (Math.Abs(maturity) < Const.QL_EPSILON)
         {
            v = 0.0;
         }
         else if (Kappa < Math.Sqrt(Const.QL_EPSILON))
         {
            v = Sigma * B(maturity, bondMaturity) * Math.Sqrt(maturity);
         }
         else
         {
            v = Sigma * B(maturity, bondMaturity) *
                Math.Sqrt(0.5 * (1.0 - Math.Exp(-2.0 * Kappa * maturity)) / Kappa);
         }
         double f = this.DiscountBond(0.0, bondMaturity, r0_);
         double k = this.DiscountBond(0.0, maturity, r0_) * strike;

         return Utils.blackFormula(type, k, f, v);
      }

      //! Short-rate dynamics in the %Vasicek model
      /*! The short-rate follows an Ornstein-Uhlenbeck process with mean
          \f$ b \f$.
      */
      public override ShortRateModel.Dynamics dynamics()
      {
         return new OneFactorModel.Dynamics(new OrnsteinUhlenbeckProcess(r0_, Kappa, Theta, Sigma));
      }
   }

   public class Gaussian : Vasicek
   {
      public Gaussian(double kappa = 0.1, double sigma = 0.01)
         : base(0, kappa, 0.0, sigma)
      {
         arguments_.Remove(arguments_[1]);
         //Ainsi, argument[2] qui est sigma devient l'argument[1]. vérifié OK!
      }
      //public override double Kappa{get { return base.Kappa; }}
      public override double Sigma
      {
         get { return arguments_[1].value(0.0); }
      }
      public override double Theta
      {
         get { return 0.0; }
      }
   }

}

