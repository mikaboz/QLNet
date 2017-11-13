using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QLNet
{
   public class BlackTimeDeterministVarianceCurve : BlackVarianceTermStructure
   {
      public VolatilityModel model_;
      public BlackTimeDeterministVarianceCurve(VolatilityModel model, Date referenceDate,Calendar calendar,DayCounter dayCounter) :
         base(referenceDate,calendar,BusinessDayConvention.Following,dayCounter)
      {
         model_ = model;
      }
      public override Date maxDate()
      {
         return Date.maxDate();
      }
      public override double maxStrike()
      {
         return double.MaxValue;
      }
      public override double minStrike()
      {
         return double.MinValue;
      }
      protected override double blackVarianceImpl(double t, double strike)
      {
         return model_.IntegratedSquareValue(0, t);
      }
   }

   /*
   public abstract class TimeDependantParameter : Parameter
   {
      public TimeDependantParameter(double[] parameters, Func<double[],double,double> function, Constraint constraint ) :
         base(parameters.Length, new Impl(function), constraint)
      {
         for (int i = 0; i < params_.Count; i++)
            params_[i] = parameters[i];
      }
      public new class Impl : Parameter.Impl
      {
         Func<double[], double, double> function_;
         public Impl(Func<double[],double,double> function)
         {
            function_ = function;
         }
         public override double value(Vector p, double t)
         {
            return function_(p.ToArray(), t);
         }
      }
   }
   */

   public class BlackScholesModel : CalibratedModel
   {
      Handle<Quote>  s0_;
      Handle<YieldTermStructure> riskFreeTermStructure_;
      Handle<YieldTermStructure> dividendeYieldTermStructure_;
      VolatilityModel volatilityModel_;
      public BlackScholesModel(Handle<Quote> s0, Handle<YieldTermStructure> riskFreeTermStructure, Handle<YieldTermStructure> dividendeYieldTermStructure, Handle<BlackTimeDeterministVarianceCurve> volatilityTS) :
         this(s0,riskFreeTermStructure,dividendeYieldTermStructure,volatilityTS.link.model_)
      {}
      public BlackScholesModel(Handle<Quote> s0, Handle<YieldTermStructure> riskFreeTermStructure, Handle<YieldTermStructure> dividendeYieldTermStructure, VolatilityModel volatilityModel) :
         base(volatilityModel.Arguments.Count)
      {
         s0_ = s0;
         volatilityModel_ = volatilityModel;
         for (int i = 0; i < arguments_.Count; i++)
            arguments_[i] = volatilityModel_.Arguments[i];

         riskFreeTermStructure_ = riskFreeTermStructure;
         dividendeYieldTermStructure_ = dividendeYieldTermStructure;

         s0_.registerWith(update);
         riskFreeTermStructure.registerWith(update);
         dividendeYieldTermStructure.registerWith(update);
      }
   }
}
