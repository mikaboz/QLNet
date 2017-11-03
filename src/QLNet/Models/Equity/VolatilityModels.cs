using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QLNet
{
   public abstract class VolatilityModel : CalibratedModel
   {
      protected VolatilityModel(int nArguments) :
         base(nArguments)
      { }
      public abstract double Value(double t);
      public abstract double IntegratedSquareValue(double t, double T);
   }
   public class ConstantVolatilityModel : VolatilityModel
   {
      public double Sigma { get { return arguments_[0].value(0.0); } }
      public ConstantVolatilityModel(double sigma) :
         base(1)
      {
         arguments_[0] = new ConstantParameter(sigma, new PositiveConstraint());
      }
      public override double IntegratedSquareValue(double t, double T)
      {
         return Sigma * Sigma * (T - t);
      }
      public override double Value(double t)
      {
         return Sigma;
      }

   }
   public abstract class TimeDeterministVolatilityModel : VolatilityModel
   {
      public TimeDeterministVolatilityModel(int nArguments) :
         base(nArguments)
      { }
   }
   public class InverseVolatilityModel : TimeDeterministVolatilityModel
   {
      private double A { get { return arguments_[0].value(0.0); } }
      private double B { get { return arguments_[1].value(0.0); } }
      private double C { get { return arguments_[2].value(0.0); } }

      public InverseVolatilityModel(double a, double b, double c) :
         base(3)
      {
         arguments_[0] = new ConstantParameter(a, new PositiveConstraint());
         arguments_[1] = new ConstantParameter(b, new PositiveConstraint());
         arguments_[2] = new ConstantParameter(c, new PositiveConstraint());
      }
      public override double Value(double t)
      {
         double t_ = (t == 0) ? (t + 0.0001) : t ;
         return A + B / Math.Pow(t_ , C);
      }
      public override double IntegratedSquareValue(double t_, double T_)
      {
         double t = (t_ == 0) ? (t_ + 0.0001) : t_;
         double T = (T_ == 0) ? (T_ + 0.0001) : T_;
         return A * (T - t) + 2 * B * A * (Math.Pow(T, 1 - C) - Math.Pow(t, 1 - C)) / (1 - C) + B * B * (Math.Pow(T, 1 - 2 * C) - Math.Pow(t, 1 - 2 * C)) / (1 - 2 * C);

      }
   }
   public class ExponentialVolatilityModel : TimeDeterministVolatilityModel
   {
      private double A { get { return arguments_[0].value(0.0); } }
      private double alpha { get { return arguments_[1].value(0.0); } }
      private double B { get { return arguments_[2].value(0.0); } }
      private double beta { get { return arguments_[3].value(0.0); } }

      public ExponentialVolatilityModel(double a, double alpha, double b, double beta) :
         base(4)
      {
         arguments_[0] = new ConstantParameter(a, new PositiveConstraint());
         arguments_[1] = new ConstantParameter(alpha, new PositiveConstraint());
         arguments_[2] = new ConstantParameter(b, new PositiveConstraint());
         arguments_[3] = new ConstantParameter(beta, new PositiveConstraint());
      }
      public override double Value(double t)
      {
         return A * Math.Exp(-alpha * t) + B * Math.Exp(-beta * t);
      }

      // Probleme: pas toujours positif ...
      public override double IntegratedSquareValue(double t, double T)
      {
         return  A*A * (Math.Exp(-2 * alpha * t) - Math.Exp(-2 * alpha * T)) / (2 * alpha)
                + 2 * A * B * (Math.Exp(-(alpha + beta) * t) - Math.Exp(-(alpha + beta) * T)) / (alpha + beta)
                + B*B * (Math.Exp(-2 * beta * t) - Math.Exp(-2 * beta * T)) / (2 * beta);
      }
   }
}
