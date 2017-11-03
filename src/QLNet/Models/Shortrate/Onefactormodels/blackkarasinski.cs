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

   public class BlackKarasinski : OneFactorModel, ITermStructureConsistentModel
   {
      public double Kappa { get { return arguments_[0].value(0.0); } }
      public double Sigma { get { return arguments_[1].value(0.0); } }

      public BlackKarasinski(Handle<YieldTermStructure> termStructure, double kappa, double sigma)
          : base(2)
      {
         arguments_[0] = new ConstantParameter(kappa, new PositiveConstraint());
         arguments_[1] = new ConstantParameter(sigma, new PositiveConstraint());
         termStructure_ = termStructure;
         termStructure.registerWith(update);
      }
      public BlackKarasinski(Handle<YieldTermStructure> termStructure)
          : this(termStructure, 0.1, 0.1)
      { }

      public override Lattice tree(TimeGrid grid)
      {
         TermStructureFittingParameter phi = new TermStructureFittingParameter(termStructure_);

         Dynamics numericDynamics =
                                 new Dynamics(phi, Kappa, Sigma);

         TrinomialTree trinomial =
                          new TrinomialTree(numericDynamics.Process, grid);
         ShortRateTree numericTree =
                          new ShortRateTree(trinomial, numericDynamics, grid);

         TermStructureFittingParameter.NumericalImpl impl =
                (TermStructureFittingParameter.NumericalImpl)phi.implementation();
         impl.reset();
         double value = 1.0;
         double vMin = -50.0;
         double vMax = 50.0;
         for (int i = 0; i < (grid.size() - 1); i++)
         {
            double discountBond = termStructure_.link.discount(grid[i + 1]);
            double xMin = trinomial.underlying(i, 0);
            double dx = trinomial.dx(i);
            Helper finder = new Helper(i, xMin, dx, discountBond, numericTree);
            Brent s1d = new Brent();
            s1d.setMaxEvaluations(1000);
            value = s1d.solve(finder, 1e-7, value, vMin, vMax);
            impl.setvalue(grid[i], value);
         }
         return numericTree;
      }

      public override ShortRateModel.Dynamics dynamics()
      {
         throw new NotImplementedException("no defined process for Black-Karasinski");
      }

      #region ITermStructureConsistentModel
      private Handle<YieldTermStructure> termStructure_;
      public Handle<YieldTermStructure> TermStructure { get { return termStructure_; } }

      #endregion

      //! Short-rate dynamics in the Black-Karasinski model
      public new class Dynamics : OneFactorModel.Dynamics
      {
         public Dynamics(Parameter fitting, double kappa, double sigma)
         : base(new OrnsteinUhlenbeckProcess(kappa, sigma))
         {
            fitting_ = fitting;
         }

         public override double variable(double t, double r)
         {
            return Math.Log(r) - fitting_.value(t);
         }

         public override double shortRate(double t, double x)
         {
            return Math.Exp(x + fitting_.value(t));
         }

         private Parameter fitting_;
      }

      // Private function used by solver to determine time-dependent parameter
      public class Helper : ISolver1d
      {
         private int size_;
         private double dt_;
         private double xMin_, dx_;
         private Vector statePrices_;
         private double discountBondPrice_;

         public Helper(int i, double xMin, double dx,
                double discountBondPrice,
                OneFactorModel.ShortRateTree tree)
         {
            size_ = tree.size(i);
            dt_ = tree.timeGrid().dt(i);
            xMin_ = xMin;
            dx_ = dx;
            statePrices_ = tree.statePrices(i);
            discountBondPrice_ = discountBondPrice;
         }

         public override double value(double theta)
         {
            double value = discountBondPrice_;
            double x = xMin_;
            for (int j = 0; j < size_; j++)
            {
               double discount = Math.Exp(-Math.Exp(theta + x) * dt_);
               value -= statePrices_[j] * discount;
               x += dx_;
            }
            return value;
         }
      }
   }
}
