/*
 Copyright (C) 2010 Philippe double (ph_real@hotmail.com)
  
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
using System.Collections.Generic;

namespace QLNet
{

   //! Two-additive-factor gaussian model class.
   /*! This class implements a two-additive-factor model defined by
       \f[
           dr_t = \varphi(t) + x_t + y_t
       \f]
       where \f$ x_t \f$ and \f$ y_t \f$ are defined by
       \f[
           dx_t = -a x_t dt + \sigma dW^1_t, x_0 = 0
       \f]
       \f[
           dy_t = -b y_t dt + \sigma dW^2_t, y_0 = 0
       \f]
       and \f$ dW^1_t dW^2_t = \rho dt \f$.

       \bug This class was not tested enough to guarantee
            its functionality.

       \ingroup shortrate
   */
   public class G2 : TwoFactorAffineModel, ITermStructureConsistentModel
   {
      #region Fitting
      TermStructureFittingParameter phi_;
      public class FittingParameter : TermStructureFittingParameter
      {
         private new class Impl : TermStructureFittingParameter.Impl
         {
            public Impl(G2 model) :
               base(model)
            { }
            public override double value(Vector v, double t)
            {
               double forward = model_.TermStructureInitialForwardRate(t);

               G2 g2 = model_ as G2;
               return forward + g2.ValueForFitting(t);
            }
         }
         public FittingParameter(G2 model) :
            base(new Impl(model))
         { }
      }
      private double ValueForFitting(double t)
      {
         double temp1 = Sigma1 * (1.0 - Math.Exp(-Kappa1 * t)) / Kappa1;
         double temp2 = Sigma2 * (1.0 - Math.Exp(-Kappa2 * t)) / Kappa2;
         return 0.5 * temp1 * temp1 + 0.5 * temp2 * temp2 +
             Rho * temp1 * temp2;
      }
      #endregion

      #region ITermStructureConsistentModel implémentation
      protected Handle<YieldTermStructure> termStructure_;
      public Handle<YieldTermStructure> TermStructure { get { return termStructure_; } }
      #endregion

      #region Accessors
      public new Gaussian First { get { return (Gaussian)base.First; } }
      public new Gaussian Second { get { return (Gaussian)base.Second; } }
      double Kappa1 { get { return First.Kappa; } }
      double Sigma1 { get { return First.Sigma; } }
      double Kappa2 { get { return Second.Kappa; } }
      double Sigma2 { get { return Second.Sigma; } }
      #endregion

      #region IAffineModel implémentation
      private double V(double t)
      {
         double expat = Math.Exp(-Kappa1 * t);
         double expbt = Math.Exp(-Kappa2 * t);
         double cx = Sigma1 / Kappa1;
         double cy = Sigma2 / Kappa2;
         double valuex = cx * cx * (t + (2.0 * expat - 0.5 * expat * expat - 1.5) / Kappa1);
         double valuey = cy * cy * (t + (2.0 * expbt - 0.5 * expbt * expbt - 1.5) / Kappa2);
         double value = 2.0 * Rho * cx * cy * (t + (expat - 1.0) / Kappa1
                                          + (expbt - 1.0) / Kappa2
                                          - (expat * expbt - 1.0) / (Kappa1 + Kappa2));
         return valuex + valuey + value;
      }
      public override double A(double t, double T)
      {
         return termStructure_.link.discount(T) / termStructure_.link.discount(t) *
          Math.Exp(0.5 * (V(T - t) - V(T) + V(t)));
      }
      public override double DiscountBondOption(Option.Type type, double strike, double maturity, double bondMaturity)
      {
         double v = sigmaP(maturity, bondMaturity);
         double f = termStructure_.link.discount(bondMaturity);
         double k = termStructure_.link.discount(maturity) * strike;
         return Utils.blackFormula(type, k, f, v);
      }
      // On préfère utiliser cette formulation pour le discount en 0, qui est moins couteuse que celle procurée par IAffineModel<TwoFactorModel>
      public override double Discount(double t)
      {
         return termStructure_.link.discount(t);
      }
      #endregion

      #region Constructors
      public G2(Handle<YieldTermStructure> termStructure, Gaussian first, Gaussian second, double rho) :
         base(first, second,rho)
      {
         termStructure_ = termStructure;
         termStructure.registerWith(update);
         generateArguments();
      }
      public G2(Handle<YieldTermStructure> termStructure, double a, double sigma, double b, double eta, double rho) :
         this(termStructure, new HullWhite(termStructure, a, sigma), new HullWhite(termStructure, b, eta), rho)
      { }
      public G2(Handle<YieldTermStructure> termStructure, double a, double sigma, double b, double eta) :
         this(termStructure, a, sigma, b, eta, -0.75)
      { }
      public G2(Handle<YieldTermStructure> termStructure, double a, double sigma, double b) :
       this(termStructure, a, sigma, b, 0.01, -0.75)
      { }
      public G2(Handle<YieldTermStructure> termStructure, double a, double sigma) :
       this(termStructure, a, sigma, 0.1, 0.01, -0.75)
      { }
      public G2(Handle<YieldTermStructure> termStructure, double a) : this(termStructure, a, 0.01, 0.1, 0.01, -0.75)
      { }
      public G2(Handle<YieldTermStructure> termStructure) :
       this(termStructure, 0.1, 0.01, 0.1, 0.01, -0.75)
      { }
      protected override void generateArguments()
      {
         phi_ = new FittingParameter(this);
      }
      #endregion

      #region Dynamics
      public override ShortRateModel.Dynamics dynamics()
      {
         return new G2.Dynamics(this);
      }
      public new class Dynamics : TwoFactorModel.Dynamics
      {
         Parameter fitting_;
         public Dynamics(G2 model) :
            base(model)
         {
            fitting_ = model.phi_;
         }
         public override double ShortRate(double t, double x, double y)
         {
            return fitting_.value(t) + x + y;
         }
      }
      #endregion

      #region Swaption valuation
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
         double start = dayCounter.yearFraction(settlement, arguments.floatingResetDates[0]);
         double w = (arguments.type == VanillaSwap.Type.Payer ? 1 : -1);
         List<double> fixedPayTimes = new InitializedList<double>(arguments.fixedPayDates.Count);
         for (int i = 0; i < fixedPayTimes.Count; ++i)
            fixedPayTimes[i] = dayCounter.yearFraction(settlement, arguments.fixedPayDates[i]);
         SwaptionPricingFunction function = new SwaptionPricingFunction(Kappa1,
                                                 Sigma1, Kappa2, Sigma2, Rho,
                                                 w, start,
                                                 fixedPayTimes,
                                                 fixedRate, this);
         double upper = function.mux() + range * function.sigmax();
         double lower = function.mux() - range * function.sigmax();
         SegmentIntegral integrator = new SegmentIntegral(intervals);
         return arguments.nominal * w * termStructure_.link.discount(start) *
             integrator.value(function.value, lower, upper);
      }
      public class SwaptionPricingFunction {
      
        #region private fields 
        double a_, sigma_, b_, eta_, rho_, w_;
        double T_;
        List<double> t_;
        double rate_;
        int size_;
        Vector A_, Ba_, Bb_;
        double mux_, muy_, sigmax_, sigmay_, rhoxy_;
        #endregion

        public SwaptionPricingFunction(double a, double sigma, double b, double eta, double rho, double w,
                                       double start, List<double> payTimes, double fixedRate, G2 model)
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

            A_  = new Vector(size_);
            Ba_ = new Vector(size_);
            Bb_ = new Vector(size_); 


            sigmax_ = sigma_*Math.Sqrt(0.5*(1.0-Math.Exp(-2.0*a_*T_))/a_);
            sigmay_ =   eta_*Math.Sqrt(0.5*(1.0-Math.Exp(-2.0*b_*T_))/b_);
            rhoxy_ = rho_*eta_*sigma_*(1.0 - Math.Exp(-(a_+b_)*T_))/
                ((a_+b_)*sigmax_*sigmay_);

            double temp = sigma_*sigma_/(a_*a_);
            mux_ = -((temp+rho_*sigma_*eta_/(a_*b_))*(1.0 - Math.Exp(-a*T_)) -
                     0.5*temp*(1.0 - Math.Exp(-2.0*a_*T_)) -
                     rho_*sigma_*eta_/(b_*(a_+b_))*
                     (1.0- Math.Exp(-(b_+a_)*T_)));

            temp = eta_*eta_/(b_*b_);
            muy_ = -((temp+rho_*sigma_*eta_/(a_*b_))*(1.0 - Math.Exp(-b*T_)) -
                     0.5*temp*(1.0 - Math.Exp(-2.0*b_*T_)) -
                     rho_*sigma_*eta_/(a_*(a_+b_))*
                     (1.0- Math.Exp(-(b_+a_)*T_)));

            for (int i = 0; i < size_; i++)
            {
               A_[i] = model.A(T_, t_[i]);
               Ba_[i] = model.First.B(T_, t_[i]); // modif: original : model.B(a_, t_[i]-T_);
               Bb_[i] = model.Second.B(T_, t_[i]); // modif: original : model.B(b_, t_[i]-T_);
            }
         }
        internal double mux() { return mux_; }
        internal double sigmax() { return sigmax_; }
        public double value(double x)  
        {
            CumulativeNormalDistribution phi = new CumulativeNormalDistribution();
            double temp = (x - mux_)/sigmax_;
            double txy = Math.Sqrt(1.0 - rhoxy_*rhoxy_);

            Vector lambda = new Vector(size_);
            int i;
            for (i=0; i<size_; i++) {
                double tau = (i==0 ? t_[0] - T_ : t_[i] - t_[i-1]);
                double c = (i==size_-1 ? (1.0+rate_*tau) : rate_*tau);
                lambda[i] = c*A_[i]*Math.Exp(-Ba_[i]*x);
            }

            SolvingFunction function = new SolvingFunction(lambda, Bb_);
            Brent s1d = new Brent();
            s1d.setMaxEvaluations(1000);
            double yb = s1d.solve(function, 1e-6, 0.00, -100.0, 100.0);

            double h1 = (yb - muy_)/(sigmay_*txy) -
                rhoxy_*(x  - mux_)/(sigmax_*txy);
            double value = phi.value(-w_*h1);


            for (i=0; i<size_; i++) {
                double h2 = h1 +
                    Bb_[i]*sigmay_*Math.Sqrt(1.0-rhoxy_*rhoxy_);
                double kappa = - Bb_[i] *
                    (muy_ - 0.5*txy*txy*sigmay_*sigmay_*Bb_[i] +
                     rhoxy_*sigmay_*(x-mux_)/sigmax_);
                value -= lambda[i] *Math.Exp(kappa)*phi.value(-w_*h2);
            }

            return Math.Exp(-0.5*temp*temp)*value/
                (sigmax_*Math.Sqrt(2.0*QLNet.Const.M_PI));
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

            public override double value(double y) {
                double value = 1.0;
                for (int i=0; i<lambda_.size(); i++) {
                    value -= lambda_[i]*Math.Exp(-Bb_[i]*y);
                }
                return value;
            }

        }
    }
      #endregion
   }
}


