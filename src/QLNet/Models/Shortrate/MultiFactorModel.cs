using System;
using System.Collections.Generic;
using System.Linq;

namespace QLNet
{
   // TODO: Tree
   public abstract class MultiFactorModel : ShortRateModel
   {
      protected bool IsCorrelatedModel;
      protected MultiFactorModel(IEnumerable<OneFactorModel> factors, double[,] correlations = null) :
         base(factors)
      {
         factors_ = factors.ToList();
         int index = 0;
         for (int i = 0; i < nFactors; i++)
            for (int k = 0; k < nArgumentsOfFactor(i); k++)
               arguments_[index++] = factors_[i].Arguments[k];
         if (correlations == null)
            IsCorrelatedModel = false;
         else
         {
            Utils.QL_REQUIRE(correlations.GetLength(0) == nFactors && correlations.GetLength(1) == nFactors, () => "matrice de correlation non adaptée en taille");
            IsCorrelatedModel = true;
            for (int i = 0; i < nFactors; i++)
               for (int j = i + 1; j < nFactors; j++)
                  arguments_.Add(new ConstantParameter(correlations[i, j], new BoundaryConstraint(-1.0, 1.0)));
         }
      }
      public int nFactors { get { return factors_.Count; } }
      private int nArgumentsOfFactor(int factorNumber) { return factors_[factorNumber].Arguments.Count; }
      public List<OneFactorModel> factors_;
      public Parameter Cor(int i, int j)
      {
         Utils.QL_REQUIRE(i < nFactors && j < nFactors, () => "OutOfRange");
         if (i == j)
            return new FixedParameter(1.0);
         else if (i > j)
            return Cor(j, i);
         else
         {
            int initialIndex = factors_.Sum(x => x.Arguments.Count) - 1; // On somme sur chaque facteur le nombre d'arguments individuels
            return arguments_[initialIndex + i + j];
         }
      }
      protected double Rho(int i, int j)
      {
         return Cor(i, j).value(0.0);
      }
      protected Matrix Cor()
      {
         if (!IsCorrelatedModel)
            return new NullCorrelationMatrix(nFactors);
         else
         {
            Matrix cor = new Matrix(nFactors, nFactors);
            for (int i = 0; i < nFactors; i++)
               for (int j = 0; j < nFactors; j++)
                  cor[i, j] = Rho(i, j);
            return cor;
         }
      }

      // TODO
      public override Lattice tree(TimeGrid t)
      {
         throw new NotImplementedException();
      }

      public override ShortRateModel.Dynamics dynamics()
      {
         return new MultiFactorModel.Dynamics(this);
      }
      public new class Dynamics : ShortRateModel.Dynamics
      {
         public List<StochasticProcess1D> Process { get; private set; }
         public Matrix Cor { get; private set; }

         public Dynamics(MultiFactorModel model)
         {
            Process = new List<StochasticProcess1D>();
            foreach (OneFactorModel mod in model.factors_)
            {
               OneFactorModel.Dynamics dynamics = (OneFactorModel.Dynamics)mod.dynamics(); //Cast pour pouvoir acceder à l'acceseur "Process"
               Process.Add(dynamics.Process);
            }
            Cor = model.Cor();
         }
         public override double shortRate(double t, Vector variables)
         {
            return variables.Sum();
         }
         /// <summary>
         /// Not Implemented for multiFactors
         /// </summary>
         public override double variable(double t, double r)
         {
            throw new NotImplementedException();
         }
      }
   }

   public abstract class MultiFactorAffineModel : MultiFactorModel, IAffineShortRateModel
   {
      protected MultiFactorAffineModel(IEnumerable<OneFactorAffineModel> factors, double[,] correlations = null) :
         base(factors, correlations)
      { }

      public abstract double A(double t, double T);
      public Vector Bvect(double t, double T)
      {
         Vector v = new Vector(nFactors);
         for (int i = 0; i < nFactors; i++)
            v[i] = ((OneFactorAffineModel)factors_[i]).B(t, T);
         return v;
      }

      public virtual double DiscountBondOption(Option.Type type, double strike, double maturity, double bondMaturity)
      {
         throw new NotImplementedException();
      }
      public double Discount(double t)
      {
         return ((IAffineShortRateModel)this).Discount(t);
      }
      public double DiscountBond(double t, double T, Vector factors)
      {
         return ((IAffineShortRateModel)this).DiscountBond(t, T, factors);
      }
   }
}
