using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;
using System.Linq;

using QLNet;

namespace TestCalibration
{
   /*
   public interface IValidatingEngine
   {
      CalibratedModel Model { get; }
      List<CalibrationHelper> Helper { get; }
      IValidatingResults getResults();
      void calculate();
   }
   public interface IValidatingResults
   {
      void reset();
   }

   // template base class for option pricing engines
   // Derived engines only need to implement the <tt>calculate()</tt> method.
   public abstract class GenericValidatingEngine<ResultsType> : IValidatingEngine
      where ResultsType : IValidatingResults, new()
   {
      protected ResultsType results_ = FastActivator<ResultsType>.Create();

      public IValidatingResults getResults()
      {
         return results_;
      }
      public void reset()
      {
         results_.reset();
      }

      public virtual void calculate()
      {
         throw new NotSupportedException();
      }
   }
   public class Reproductibilite : GenericValidatingEngine<Reproductibilite.Result>
   {
      CalibratedModel model_;
      List<CalibrationHelper> helpers_;
      public CalibratedModel Model { get { return model_; } }
      public List<CalibrationHelper> Helpers { get { return helpers_; } }
      public Reproductibilite(CalibratedModel model, List<CalibrationHelper> helpers)
      {
         model_ = model;
         helpers_ = helpers;
      }

      public override void calculate()
      {
         foreach (CalibrationHelper helper in helpers_)
         {
            helper.SetCalibrationErrorType(CalibrationHelper.CalibrationErrorType.RelativePriceError);
            double calibrationError = helper.calibrationError();
            double individualSSE = calibrationError * calibrationError;
            results_.SSE_.Add(individualSSE);
         }
      }

      public class Result : IValidatingResults
      {
         public List<double> SSE_;
         public Result()
         {
            SSE_ = new List<double>();
         }
         public void reset()
         {
            SSE_.Clear();
         }
      }

      public static class ExtensionToStringIEnumerable
      {
         public static string ToString<T>(this IEnumerable<T> enumerable)
         {
            string str = "";
            foreach (T t in enumerable)
               str += t.ToString() + "\t";
            return str;
         }
      }

      [TestClass]
      public class T_BlackScholesCalibration
      {
         public class Validation
         {
            private double tolerance_;
            private CalibrationHelper.CalibrationErrorType errorType_;
            public Validation(double tolerance, CalibrationHelper.CalibrationErrorType errorType)
            {
               errorType_ = errorType;
               tolerance_ = tolerance;
            }
            public void Validate(CalibratedModel model, List<CalibrationHelper> helpers)
            {
               List<double> SSE = new List<double>();
               foreach (CalibrationHelper helper in helpers)
               {
                  helper.SetCalibrationErrorType(errorType_);
                  List<CalibrationHelper> individual = new List<CalibrationHelper>();
                  individual.Add(helper);
                  double calibrationError = helper.calibrationError();
                  double individualSSE = calibrationError * calibrationError;
                  SSE.Add(individualSSE);
               }

               double CalibrationFunctionValue = Math.Sqrt(SSE.Sum());
               bool test = CalibrationFunctionValue < tolerance_;

               Assert.IsTrue(test,
                  "\n individualRelativeErrors: " + SSE.ToString<double>() +
                  "\ntolerance: " + tolerance_ +
                  "\nvalue of calibrationFunction on helpers with uniform weights: " + CalibrationFunctionValue +
                  "\nLe modèle n'a pas réussi à reproduire les prix de marché conformément à la tolérance exigée"
                  );
            }
         }
         [TestMethod]
         public void ConstantVSConstant()
         {
            #region Global
            Calendar calendar = new TARGET();
            Date referenceDate = DateTime.Now;
            DayCounter dayCounter = new Actual365Fixed();
            #endregion

            #region Stock
            double spot_ = 100;
            double flatDividend_ = 0.0;
            Handle<Quote> spot = new Handle<Quote>(new SimpleQuote(spot_));
            Handle<YieldTermStructure> dividendTS = new Handle<YieldTermStructure>(new FlatForward(referenceDate, flatDividend_, dayCounter));
            #endregion

            #region RiskFree
            double flatRiskFreeRate_ = 0.01;
            Handle<YieldTermStructure> riskFreeTS = new Handle<YieldTermStructure>(new FlatForward(referenceDate, flatRiskFreeRate_, dayCounter));
            #endregion

            #region CalibrationHelpers
            int[] timeToMaturity = new int[10];
            for (int i = 0; i < timeToMaturity.Length; i++)
            {
               timeToMaturity[i] = i + 1;
            }
            TimeUnit timeUnit = TimeUnit.Years;
            double[] strike = new double[15];
            for (int i = 0; i < strike.Length; i++)
            {
               strike[i] = 50 + 10 * i;
            }
            double constantVol = 0.1;
            //Console.WriteLine("generatorsParameters:\nconstantVol: {0}", constantVol);
            Handle<BlackVolTermStructure> blackConstantVol = new Handle<BlackVolTermStructure>(new BlackConstantVol(referenceDate, calendar, constantVol, dayCounter));
            BlackScholesMertonProcess generatorBsmProcess = new BlackScholesMertonProcess(spot, dividendTS, riskFreeTS, new Handle<BlackVolTermStructure>(blackConstantVol));
            IPricingEngine generatorEngine = new AnalyticEuropeanEngine(generatorBsmProcess);
            List<CalibrationHelper> helpers = new List<CalibrationHelper>();
            for (int i = 0; i < timeToMaturity.Length; i++)
               for (int j = 0; j < strike.Length; j++)
               {
                  Period maturity = new Period(timeToMaturity[i], timeUnit);
                  StrikedTypePayoff payOff = new PlainVanillaPayoff(Option.Type.Call, strike[j]);
                  Exercise exercise = new EuropeanExercise(Settings.evaluationDate() + maturity);
                  VanillaOption option = new VanillaOption(payOff, exercise);
                  option.setPricingEngine(generatorEngine);
                  double price = option.NPV();
                  double impliedVol = option.impliedVolatility(price, generatorBsmProcess);
                  Handle<Quote> vol = new Handle<Quote>(new SimpleQuote(impliedVol));
                  CalibrationHelper helper = new VanillaOptionHelper(maturity, calendar, spot, strike[j], vol, riskFreeTS, dividendTS);
                  helpers.Add(helper);
               }
            #endregion

            #region Model

            #region Volatility
            double guessConstantVolatility = 0.2;
            Console.WriteLine("Seed:\nguessConstantVolatility: {0}", guessConstantVolatility);
            IParametricVolatility guessConstantVolatilityParameter = new ConstantVolatilityParameter(guessConstantVolatility);
            ParametricVolatilityTermStructure volatilityTS = new ParametricVolatilityTermStructure(guessConstantVolatilityParameter, referenceDate, dayCounter, calendar);
            #endregion

            BlackScholesModel model = new BlackScholesModel(spot, dividendTS, riskFreeTS, volatilityTS);
            GeneralizedBlackScholesProcess bsmProcess = model.process_;
            IPricingEngine engine = new AnalyticEuropeanEngine(bsmProcess);
            foreach (CalibrationHelper c in helpers)
               c.setPricingEngine(engine);
            #endregion

            #region Optimizer
            OptimizationMethod LM = new LevenbergMarquardt();
            EndCriteria endCriteria = new EndCriteria(100, 5, 10e-4, 10e-4, 10e-5);
            #endregion

            #region Calibration
            model.calibrate(helpers, LM, endCriteria);
            #endregion

            #region Results
            //Console.WriteLine("output");
            double[] output = model.parameters().ToArray();
            //foreach (double d in output)
            // Console.WriteLine(d);
            //Console.ReadLine();
            #endregion

            #region TEST
            Validation validation = new Validation(0.02, CalibrationHelper.CalibrationErrorType.RelativePriceError);
            validation.Validate(model, helpers);
            #endregion
         }

         [TestMethod]
         public void ConstantVSExponential()
         {
            #region Global
            Calendar calendar = new TARGET();
            Date referenceDate = DateTime.Now;
            DayCounter dayCounter = new Actual365Fixed();
            #endregion

            #region Stock
            double spot_ = 100;
            double flatDividend_ = 0.0;
            Handle<Quote> spot = new Handle<Quote>(new SimpleQuote(spot_));
            Handle<YieldTermStructure> dividendTS = new Handle<YieldTermStructure>(new FlatForward(referenceDate, flatDividend_, dayCounter));
            #endregion

            #region RiskFree
            double flatRiskFreeRate_ = 0.01;
            Handle<YieldTermStructure> riskFreeTS = new Handle<YieldTermStructure>(new FlatForward(referenceDate, flatRiskFreeRate_, dayCounter));
            #endregion

            #region CalibrationHelpers

            int[] timeToMaturity = new int[10];
            for (int i = 0; i < timeToMaturity.Length; i++)
            {
               timeToMaturity[i] = i + 1;
            }
            TimeUnit timeUnit = TimeUnit.Years;
            double[] strike = new double[15];
            for (int i = 0; i < strike.Length; i++)
            {
               strike[i] = 50 + 10 * i;
            }
            double A = 0.2;
            double alpha = 0.1;
            double B = 0.9;
            double beta = 0.2;
            Period horizon = new Period(timeToMaturity.Max(), timeUnit);
            IParametricVolatility exponentialVolatilityParameter = new ExponentialVolatilityParameter(A, alpha, B, beta, horizon);
            ParametricVolatilityTermStructure exponentialVolatilityTS = new ParametricVolatilityTermStructure(exponentialVolatilityParameter, referenceDate, dayCounter, calendar);
            BlackScholesModel generatorModel = new BlackScholesModel(spot, dividendTS, riskFreeTS, exponentialVolatilityTS);
            GeneralizedBlackScholesProcess generatorBsmProcess = generatorModel.process_;
            IPricingEngine generatorEngine = new AnalyticEuropeanEngine(generatorBsmProcess);
            List<CalibrationHelper> helpers = new List<CalibrationHelper>();
            for (int i = 0; i < timeToMaturity.Length; i++)
               for (int j = 0; j < strike.Length; j++)
               {
                  Period maturity = new Period(timeToMaturity[i], timeUnit);
                  StrikedTypePayoff payOff = new PlainVanillaPayoff(Option.Type.Call, strike[j]);
                  Exercise exercise = new EuropeanExercise(Settings.evaluationDate() + maturity);
                  VanillaOption option = new VanillaOption(payOff, exercise);
                  option.setPricingEngine(generatorEngine);
                  double price = option.NPV();
                  double impliedVol = option.impliedVolatility(price, generatorBsmProcess);
                  helpers.Add(new VanillaOptionHelper(maturity, calendar, spot, strike[j], new Handle<Quote>(new SimpleQuote(impliedVol)), riskFreeTS, dividendTS));
               }
            #endregion

            #region Model

            #region Volatility
            double guessVolatility = 0.2;
            IParametricVolatility constantVolatility = new ConstantVolatilityParameter(guessVolatility);
            ParametricVolatilityTermStructure volatilityTS = new ParametricVolatilityTermStructure(constantVolatility, referenceDate, dayCounter, calendar);
            #endregion
            BlackScholesModel model = new BlackScholesModel(spot, dividendTS, riskFreeTS, volatilityTS);
            GeneralizedBlackScholesProcess bsmProcess = model.process_;
            IPricingEngine engine = new AnalyticEuropeanEngine(bsmProcess);
            foreach (CalibrationHelper c in helpers)
               c.setPricingEngine(engine);
            #endregion

            #region Optimizer
            OptimizationMethod LM = new LevenbergMarquardt();
            EndCriteria endCriteria = new EndCriteria(100, 5, 10e-4, 10e-4, 10e-5);
            #endregion

            #region Calibration
            model.calibrate(helpers, LM, endCriteria);
            #endregion

            #region Results
            double[] output = model.parameters().ToArray();
            foreach (double d in output)
               Console.WriteLine(d);
            Console.ReadLine();
            #endregion

            #region BackTest
            foreach (CalibrationHelper helper in helpers)
            {
               Console.WriteLine("marketValue: {0}\tmodelValue: {1}", helper.marketValue(), helper.modelValue());
            }
            Console.ReadLine();
            #endregion
         }
         public static void ExponentialVSConstant()
         {
            #region Global
            Calendar calendar = new TARGET();
            Date referenceDate = DateTime.Now;
            DayCounter dayCounter = new Actual365Fixed();
            #endregion

            #region Stock
            double spot_ = 100;
            double flatDividend_ = 0.0;
            Handle<Quote> spot = new Handle<Quote>(new SimpleQuote(spot_));
            Handle<YieldTermStructure> dividendTS = new Handle<YieldTermStructure>(new FlatForward(referenceDate, flatDividend_, dayCounter));
            #endregion

            #region RiskFree
            double flatRiskFreeRate_ = 0.01;
            Handle<YieldTermStructure> riskFreeTS = new Handle<YieldTermStructure>(new FlatForward(referenceDate, flatRiskFreeRate_, dayCounter));
            #endregion

            #region CalibrationHelpers

            int[] timeToMaturity = new int[10];
            for (int i = 0; i < timeToMaturity.Length; i++)
            {
               timeToMaturity[i] = i + 1;
            }
            TimeUnit timeUnit = TimeUnit.Years;
            double[] strike = new double[15];
            for (int i = 0; i < strike.Length; i++)
            {
               strike[i] = 50 + 10 * i;
            }
            double constantVol = 0.1;
            Period horizon = new Period(timeToMaturity.Max(), timeUnit);
            Handle<BlackVolTermStructure> constantVolTS = new Handle<BlackVolTermStructure>(new BlackConstantVol(referenceDate, calendar, constantVol, dayCounter));
            BlackScholesMertonProcess generatorBsmProcess = new BlackScholesMertonProcess(spot, dividendTS, riskFreeTS, constantVolTS);
            IPricingEngine generatorEngine = new AnalyticEuropeanEngine(generatorBsmProcess);
            List<CalibrationHelper> helpers = new List<CalibrationHelper>();
            for (int i = 0; i < timeToMaturity.Length; i++)
               for (int j = 0; j < strike.Length; j++)
               {
                  Period maturity = new Period(timeToMaturity[i], timeUnit);
                  StrikedTypePayoff payOff = new PlainVanillaPayoff(Option.Type.Call, strike[j]);
                  Exercise exercise = new EuropeanExercise(Settings.evaluationDate() + maturity);
                  VanillaOption option = new VanillaOption(payOff, exercise);
                  option.setPricingEngine(generatorEngine);
                  double price = option.NPV();
                  double impliedVol = option.impliedVolatility(price, generatorBsmProcess);
                  helpers.Add(new VanillaOptionHelper(maturity, calendar, spot, strike[j], new Handle<Quote>(new SimpleQuote(impliedVol)), riskFreeTS, dividendTS));
               }
            #endregion

            #region Model

            #region Volatility
            double A = 0.2;
            double alpha = 0.1;
            double B = 0.9;
            double beta = 0.2;
            IParametricVolatility exponentialVolatilityParameter = new ExponentialVolatilityParameter(A, alpha, B, beta, horizon);
            ParametricVolatilityTermStructure exponentialVolatilityTS = new ParametricVolatilityTermStructure(exponentialVolatilityParameter, referenceDate, dayCounter, calendar);
            #endregion

            BlackScholesModel model = new BlackScholesModel(spot, dividendTS, riskFreeTS, exponentialVolatilityTS);
            GeneralizedBlackScholesProcess bsmProcess = model.process_;
            IPricingEngine engine = new AnalyticEuropeanEngine(bsmProcess);
            foreach (CalibrationHelper c in helpers)
               c.setPricingEngine(engine);
            #endregion

            #region Optimizer
            OptimizationMethod LM = new LevenbergMarquardt();
            EndCriteria endCriteria = new EndCriteria(100, 30, 10e-4, 10e-4, 10e-5);
            #endregion

            #region Calibration
            model.calibrate(helpers, LM, endCriteria);
            #endregion

            #region Results
            double[] output = model.parameters().ToArray();
            foreach (double d in output)
               Console.WriteLine(d);
            Console.ReadLine();
            #endregion

            #region BackTest
            foreach (CalibrationHelper helper in helpers)
            {
               Console.WriteLine("marketValue: {0}\tmodelValue: {1}", helper.marketValue(), helper.modelValue());
            }
            Console.ReadLine();
            #endregion
         }
         public static void ExponentialVSExponential()
         {
            #region Global
            Calendar calendar = new TARGET();
            Date referenceDate = DateTime.Now;
            DayCounter dayCounter = new Actual365Fixed();
            #endregion

            #region Stock
            double spot_ = 100;
            double flatDividend_ = 0.0;
            Handle<Quote> spot = new Handle<Quote>(new SimpleQuote(spot_));
            Handle<YieldTermStructure> dividendTS = new Handle<YieldTermStructure>(new FlatForward(referenceDate, flatDividend_, dayCounter));
            #endregion

            #region RiskFree
            double flatRiskFreeRate_ = 0.01;
            Handle<YieldTermStructure> riskFreeTS = new Handle<YieldTermStructure>(new FlatForward(referenceDate, flatRiskFreeRate_, dayCounter));
            #endregion

            #region CalibrationHelpers

            int[] timeToMaturity = new int[10];
            for (int i = 0; i < timeToMaturity.Length; i++)
            {
               timeToMaturity[i] = i + 1;
            }
            TimeUnit timeUnit = TimeUnit.Years;
            double[] strike = new double[15];
            for (int i = 0; i < strike.Length; i++)
            {
               strike[i] = 50 + 10 * i;
            }
            double A = 0.2;
            double alpha = 0.6;
            double B = 0.1;
            double beta = 0.5;
            Console.WriteLine("generatorsParameters:\nA: {0}\nalpha: {1}\nB: {2}\nbeta: {3}", A, alpha, B, beta);
            Period horizon = new Period(timeToMaturity.Max(), timeUnit);
            TimeDeterministParameter generatorExpVolatilityParameter = new ExponentialVolatilityParameter(A, alpha, B, beta, horizon);
            Handle<BlackVolTermStructure> generatorExpVolatilityTS = new Handle<BlackVolTermStructure>(new ParametricVolatilityTermStructure(generatorExpVolatilityParameter, referenceDate, dayCounter, calendar));
            BlackScholesMertonProcess generatorBsmProcess = new BlackScholesMertonProcess(spot, dividendTS, riskFreeTS, generatorExpVolatilityTS);
            IPricingEngine generatorEngine = new AnalyticEuropeanEngine(generatorBsmProcess);
            List<CalibrationHelper> helpers = new List<CalibrationHelper>();
            for (int i = 0; i < timeToMaturity.Length; i++)
               for (int j = 0; j < strike.Length; j++)
               {
                  Period maturity = new Period(timeToMaturity[i], timeUnit);
                  StrikedTypePayoff payOff = new PlainVanillaPayoff(Option.Type.Call, strike[j]);
                  Exercise exercise = new EuropeanExercise(Settings.evaluationDate() + maturity);
                  VanillaOption option = new VanillaOption(payOff, exercise);
                  option.setPricingEngine(generatorEngine);
                  double price = option.NPV();
                  double impliedVol = option.impliedVolatility(price, generatorBsmProcess);
                  helpers.Add(new VanillaOptionHelper(maturity, calendar, spot, strike[j], new Handle<Quote>(new SimpleQuote(impliedVol)), riskFreeTS, dividendTS));
               }
            #endregion

            #region Model

            #region Volatility
            double guessA = 0.05;
            double guessalpha = 0.3;
            double guessB = 0.22;
            double guessbeta = 0.12;
            Console.WriteLine("seed:\nA: {0}\n: {1}\nB: {2}\nbeta: {3}", guessA, guessalpha, guessB, guessbeta);
            IParametricVolatility exponentialVolatilityParameter = new ExponentialVolatilityParameter(guessA, guessalpha, guessB, guessbeta, horizon);
            ParametricVolatilityTermStructure exponentialVolatilityTS = new ParametricVolatilityTermStructure(exponentialVolatilityParameter, referenceDate, dayCounter, calendar);
            #endregion

            BlackScholesModel model = new BlackScholesModel(spot, dividendTS, riskFreeTS, exponentialVolatilityTS);
            GeneralizedBlackScholesProcess bsmProcess = model.process_;
            Console.ReadLine();
            IPricingEngine engine = new AnalyticEuropeanEngine(bsmProcess);
            foreach (CalibrationHelper c in helpers)
               c.setPricingEngine(engine);
            #endregion

            #region Optimizer
            OptimizationMethod LM = new LevenbergMarquardt();
            EndCriteria endCriteria = new EndCriteria(100, null, 10e-4, 10e-4, null);
            #endregion

            #region Calibration
            model.calibrate(helpers, LM, endCriteria);
            #endregion

            #region Results
            double[] output = model.parameters().ToArray();
            Console.WriteLine("output: ");
            foreach (double d in output)
               Console.WriteLine(d);
            Console.ReadLine();
            #endregion

            #region BackTest
            foreach (CalibrationHelper helper in helpers)
            {
               Console.WriteLine("marketValue: {0}\tmodelValue: {1}", helper.marketValue(), helper.modelValue());
            }
            Console.ReadLine();
            #endregion
         }
      }
   }
   */
}
