using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ESG
{
   using QLNet;

   class Program
   {
      static void Main(string[] args)
      {
         BlackScholesCalibration.ExponentialVSExponential();
      }
   }

   public static class BlackScholesCalibration
   {
      public static void ConstantVSConstant()
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
         Console.WriteLine("generatorsParameters:\nconstantVol: {0}", constantVol);
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
               helpers.Add(new VanillaOptionHelper(maturity, calendar, spot, strike[j], new Handle<Quote>(new SimpleQuote(impliedVol)), riskFreeTS, dividendTS));
            }
         #endregion

         #region Model

         #region Volatility
         double guessConstantVolatility = 0.2;
         Console.WriteLine("Seed:\nguessConstantVolatility: {0}", guessConstantVolatility);
         Parameter guessConstantVolatilityParameter = new ConstantVolatilityParameter(guessConstantVolatility);
         ParametricVolatilityTermStructure volatilityTS = new ParametricVolatilityTermStructure(guessConstantVolatilityParameter, referenceDate, dayCounter, calendar);
         #endregion

         BlackScholesModel model = new BlackScholesModel(spot, dividendTS, riskFreeTS, volatilityTS);
         BlackScholesMertonProcess bsmProcess = model.process_;
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
         Console.WriteLine("output");
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
      public static void ConstantVSExponential()
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
         Parameter exponentialVolatilityParameter = new ExponentialVolatilityParameter(A, alpha, B, beta, horizon);
         ParametricVolatilityTermStructure exponentialVolatilityTS = new ParametricVolatilityTermStructure(exponentialVolatilityParameter, referenceDate, dayCounter, calendar);
         BlackScholesModel generatorModel = new BlackScholesModel(spot, dividendTS, riskFreeTS, exponentialVolatilityTS);
         BlackScholesMertonProcess generatorBsmProcess = generatorModel.process_;
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
         Parameter constantVolatility = new ConstantVolatilityParameter(guessVolatility);
         ParametricVolatilityTermStructure volatilityTS = new ParametricVolatilityTermStructure(constantVolatility, referenceDate, dayCounter, calendar);
         #endregion
         BlackScholesModel model = new BlackScholesModel(spot, dividendTS, riskFreeTS, volatilityTS);
         BlackScholesMertonProcess bsmProcess = model.process_;
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
         Parameter exponentialVolatilityParameter = new ExponentialVolatilityParameter(A, alpha, B, beta, horizon);
         ParametricVolatilityTermStructure exponentialVolatilityTS = new ParametricVolatilityTermStructure(exponentialVolatilityParameter, referenceDate, dayCounter, calendar);
         #endregion

         BlackScholesModel model = new BlackScholesModel(spot, dividendTS, riskFreeTS, exponentialVolatilityTS);
         BlackScholesMertonProcess bsmProcess = model.process_;
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
         Parameter exponentialVolatilityParameter = new ExponentialVolatilityParameter(guessA, guessalpha, guessB, guessbeta, horizon);
         ParametricVolatilityTermStructure exponentialVolatilityTS = new ParametricVolatilityTermStructure(exponentialVolatilityParameter, referenceDate, dayCounter, calendar);
         #endregion

         BlackScholesModel model = new BlackScholesModel(spot, dividendTS, riskFreeTS, exponentialVolatilityTS);
         BlackScholesMertonProcess bsmProcess = model.process_;
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

   public static class HestonCalibration
   {
      public static void HestonVSConstant()
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
         double v0 = 0.2;
         double kappa = 0.1;
         double theta = 0.9;
         double sigma = 0.2;
         double rho = 0.6;
         HestonProcess hestonProcess = new HestonProcess(riskFreeTS, dividendTS, spot, v0, kappa, theta, sigma, rho);
         HestonModel model = new HestonModel(hestonProcess);
         IPricingEngine engine = new AnalyticHestonEngine(model);
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
      public static void HestonVSHeston()
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
         double v0 = 0.4;
         double kappa = 0.2;
         double theta = 0.3;
         double sigma = 0.5;
         double rho = 0.7;
         Console.WriteLine("generatorsParameters:\nv0: {0}\nkappa: {1}\ntheta: {2}\nsigma: {3}\nrho: {4}", v0, kappa, theta, sigma, rho);
         HestonProcess generatorHestonProcess = new HestonProcess(riskFreeTS, dividendTS, spot, v0, kappa, theta, sigma, rho);
         HestonModel generatorHestonModel = new HestonModel(generatorHestonProcess);
         IPricingEngine generatorEngine = new AnalyticHestonEngine(generatorHestonModel);
         Handle<BlackVolTermStructure> hestonVolSurface = new Handle<BlackVolTermStructure>(new HestonBlackVolSurface(new Handle<HestonModel>(generatorHestonModel)));
         GeneralizedBlackScholesProcess gbsProcess = new GeneralizedBlackScholesProcess(spot, dividendTS, riskFreeTS, hestonVolSurface);
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
               double impliedVol = option.impliedVolatility(price, gbsProcess);
               helpers.Add(new VanillaOptionHelper(maturity, calendar, spot, strike[j], new Handle<Quote>(new SimpleQuote(impliedVol)), riskFreeTS, dividendTS));
            }
         #endregion

         #region Model
         double guessv0 = 0.2;
         double guesskappa = 0.1;
         double guesstheta = 0.9;
         double guesssigma = 0.2;
         double guessrho = 0.6;
         Console.WriteLine("guess:\nv0: {0}\nkappa: {1}\ntheta: {2}\nsigma: {3}\nrho: {4}", guessv0, guesskappa, guesstheta, guesssigma, guessrho);
         HestonProcess hestonProcess = new HestonProcess(riskFreeTS, dividendTS, spot, guessv0, guesskappa, guesstheta, guesssigma, guessrho);
         HestonModel model = new HestonModel(hestonProcess);
         IPricingEngine engine = new AnalyticHestonEngine(model);
         foreach (CalibrationHelper c in helpers)
            c.setPricingEngine(engine);
         #endregion

         #region Optimizer
         OptimizationMethod LM = new LevenbergMarquardt();
         EndCriteria endCriteria = new EndCriteria(100, 30, 10e-4, 10e-4, 10e-5);
         #endregion

         #region Calibration
         PenalizationFunction penalization = new HestonModel.Penalization(model, 10e8, 1);
         model.calibrate(helpers, LM, endCriteria,null,null,null,penalization);
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
   }
}
