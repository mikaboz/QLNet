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
   /// Single-factor Hull-White (extended %Vasicek) model class.
   /// <remarks>
   /// This class implements the standard single-factor Hull-White model defined by
   /// dr_t = (\theta(t) - \alpha r_t)dt + \sigma dW_t
   /// where alpha and sigma are constants.
   /// calibration results are tested against cached values
   /// When the term structure is relinked, the r0 parameter of
   /// the underlying Vasicek model is not updated.
   /// </remarks>
   /// </summary>
   public class HullWhite : Gaussian, ITermStructureConsistentModel
   {
      #region ITermStructureConsistentModel
      private Handle<YieldTermStructure> termStructure_;
      public Handle<YieldTermStructure> TermStructure { get { return termStructure_; } }
      #endregion

      private TermStructureFittingParameter phi_;
      /// <summary>
      /// r(t) = x(t) + phi(t)
      /// where dx(t) = -kappa*x(t)*dt + sigma*dW(t)
      /// </summary>
      /// <param name="termStructure">term structure to fit phi(t) with</param>
      /// <param name="kappa">time-constant coefficient </param>
      /// <param name="sigma">time-constant coefficient </param>
      public HullWhite(Handle<YieldTermStructure> termStructure, double kappa = 0.1, double sigma = 0.01)
         : base(kappa, sigma)
      {
         this.termStructure_ = termStructure;
         termStructure.registerWith(update);
         generateArguments();
      }

      public override Lattice tree(TimeGrid grid)
      {
         TermStructureFittingParameter phi = new TermStructureFittingParameter(termStructure_);
         Dynamics numericDynamics = new Dynamics(phi, Kappa, Sigma);
         TrinomialTree trinomial = new TrinomialTree( numericDynamics.Process, grid);
         ShortRateTree numericTree = new ShortRateTree(trinomial, numericDynamics, grid);
         TermStructureFittingParameter.NumericalImpl impl =
            (TermStructureFittingParameter.NumericalImpl) phi.implementation();
         impl.reset();
         for (int i = 0; i < (grid.size() - 1); i++)
         {
            double discountBond = termStructure_.link.discount(grid[i + 1]);
            Vector statePrices = numericTree.statePrices(i);
            int size = numericTree.size(i);
            double dt = numericTree.timeGrid().dt(i);
            double dx = trinomial.dx(i);
            double x = trinomial.underlying(i, 0);
            double value = 0.0;
            for (int j = 0; j < size; j++)
            {
               value += statePrices[j] * Math.Exp(-x * dt);
               x += dx;
            }
            value = Math.Log(value / discountBond) / dt;
            impl.setvalue(grid[i], value);
         }
         return numericTree;
      }
      public override ShortRateModel.Dynamics dynamics()
      {
         return new HullWhite.Dynamics(phi_, Kappa, Sigma);
      }
      public override double DiscountBondOption(Option.Type type, double strike, double maturity, double bondMaturity)
      {
         double _a = Kappa;
         double v;
         if (_a < Math.Sqrt(Const.QL_EPSILON))
         {
            v = Sigma * B(maturity, bondMaturity) * Math.Sqrt(maturity);
         }
         else
         {
            v = Sigma * B(maturity, bondMaturity) *
                Math.Sqrt(0.5 * (1.0 - Math.Exp(-2.0 * _a * maturity)) / _a);
         }
         double f = termStructure_.link.discount(bondMaturity);
         double k = termStructure_.link.discount(maturity) * strike;

         return Utils.blackFormula(type, k, f, v);
      }

      /*! Futures convexity bias (i.e., the difference between
          futures implied rate and forward rate) calculated as in
          G. Kirikos, D. Novak, "Convexity Conundrums", Risk
          Magazine, March 1997.

          \note t and T should be expressed in yearfraction using
                deposit day counter, F_quoted is futures' market price.
      */
      public static double convexityBias(double futuresPrice, double t, double T, double sigma, double a)
      {
         Utils.QL_REQUIRE(futuresPrice >= 0.0, () => "negative futures price (" + futuresPrice + ") not allowed");
         Utils.QL_REQUIRE(t >= 0.0, () => "negative t (" + t + ") not allowed");
         Utils.QL_REQUIRE(T >= t, () => "T (" + T + ") must not be less than t (" + t + ")");
         Utils.QL_REQUIRE(sigma >= 0.0, () => "negative sigma (" + sigma + ") not allowed");
         Utils.QL_REQUIRE(a >= 0.0, () => "negative a (" + a + ") not allowed");

         double deltaT = (T - t);
         double tempDeltaT = (1.0 - Math.Exp(-a * deltaT)) / a;
         double halfSigmaSquare = sigma * sigma / 2.0;

         // lambda adjusts for the fact that the underlying is an interest rate
         double lambda = halfSigmaSquare * (1.0 - Math.Exp(-2.0 * a * t)) / a *
                         tempDeltaT * tempDeltaT;

         double tempT = (1.0 - Math.Exp(-a * t)) / a;

         // phi is the MtM adjustment
         double phi = halfSigmaSquare * tempDeltaT * tempT * tempT;

         // the adjustment
         double z = lambda + phi;

         double futureRate = (100.0 - futuresPrice) / 100.0;
         return (1.0 - Math.Exp(-z)) * (futureRate + 1.0 / (T - t));
      }
      
      protected override void generateArguments()
      {
         phi_ = new FittingParameter(this);
      }

      public override double A(double t, double T)
      {
         double PM0t= termStructure_.link.discount(t);
         double PM0T = termStructure_.link.discount(T);
         double V0T = base.V(0, T);
         double V0t = base.V(0, t);
         double added = PM0T / PM0t * Math.Exp(0.5 * (V0t - V0T));
         return base.A(t, T) * added;
         /*
         double discount1 = termStructure_.link.discount(t);
         double discount2 = termStructure_.link.discount(T);

         double BtT = B(t, T);
         double temp = Sigma * BtT;
         double value = BtT * this.TermStructureInitialForwardRate(t) - 0.25 * temp * temp * B(0.0, 2.0 * t);
         return Math.Exp(value) * discount2 / discount1;
         */
      }

      //! Short-rate dynamics in the Hull-White model
      public new class Dynamics : OneFactorModel.Dynamics
      {
         private TermStructureFittingParameter fitting_;
         public Dynamics(TermStructureFittingParameter fitting, double kappa, double sigma)
            : base(new OrnsteinUhlenbeckProcess(kappa, sigma))
         {
            fitting_ = fitting;
         }
         public override double variable(double t, double r)
         {
            return r - fitting_.value(t);
         }
         public override double shortRate(double t, double x)
         {
            return x + fitting_.value(t);
         }
      
      }

      //! Analytical term-structure fitting parameter \f$ \varphi(t) \f$.
      /*! \f$ \varphi(t) \f$ is analytically defined by
          \f[
              \varphi(t) = f(t) + \frac{1}{2}[\frac{\sigma(1-e^{-at})}{a}]^2,
          \f]
          where \f$ f(t) \f$ is the instantaneous forward rate at \f$ t \f$.
      */
      private double ValueForFitting(double t)
      {
         double temp = Kappa < Math.Sqrt(Const.QL_EPSILON) ? Sigma * t : Sigma * (1.0 - Math.Exp(-Kappa * t)) / Kappa;
         return 0.5 * temp * temp;
      }
      public class FittingParameter : TermStructureFittingParameter
      {
         private new class Impl : TermStructureFittingParameter.Impl
         {
            public Impl(HullWhite model) :
               base(model)
            {}
            // Cas spécial pour HullWhite, expression de alpha différente que pour les models affines.
            public override double value(Vector v, double t)
            {
               return model_.TermStructureInitialForwardRate(t) + ((HullWhite)model_).ValueForFitting(t);
            }
         }
         public FittingParameter(HullWhite model)
            : base(new FittingParameter.Impl(model))
         { }
      }
   }
}
