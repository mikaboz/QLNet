using System;

namespace QLNet
{
   public class CIR2PlusPlus : CIR2, ITermStructureConsistentModel
   {
      #region Fitting
      public class FittingParameter : TermStructureFittingParameter
      {
         public FittingParameter(CIR2PlusPlus model) :
            base(new Impl(model))
         { }
         private new class Impl : TermStructureFittingParameter.Impl
         {
            public Impl(CIR2PlusPlus model) :
               base(model)
            { }
            public override double value(Vector p, double t)
            {
               CIR2 model = (CIR2)model_;
               CIRPlusPlus cirPP1 = (CIRPlusPlus)model.First;
               CIRPlusPlus cirPP2 = (CIRPlusPlus)model.Second;
               return model_.TermStructureInitialForwardRate(t) - cirPP1.ValueForFitting(t) - cirPP2.ValueForFitting(t);
            }
         }
      }
      private FittingParameter phi_;
      #endregion

      #region ItermStructureConsistentModel implementation
      private Handle<YieldTermStructure> termStructure_;
      public Handle<YieldTermStructure> TermStructure { get { return termStructure_; } }
      #endregion

      #region Constructors
      public CIR2PlusPlus(Handle<YieldTermStructure> termStructure, double x0, double kappa1, double theta1, double sigma1, double y0, double kappa2, double theta2, double sigma2) :
         this(termStructure, new CoxIngersollRoss(x0,kappa1, theta1, sigma1), new CoxIngersollRoss(y0, kappa2, theta2, sigma2))
      {}
      public CIR2PlusPlus(Handle<YieldTermStructure> termStructure, CoxIngersollRoss first, CoxIngersollRoss second):
         base(first, second)
      {
            termStructure_ = termStructure;
            termStructure.registerWith(update);
            generateArguments();
      }
      protected override void generateArguments()
      {
         phi_ = new FittingParameter(this);
      }
      #endregion

      #region Dynamics
      public override ShortRateModel.Dynamics dynamics()
      {
         return new CIR2PlusPlus.Dynamics(this);
      }
      public new class Dynamics : CIR2.Dynamics
      {
         private FittingParameter fitting_;
         public Dynamics(CIR2PlusPlus model) :
            base(model)
         {
            fitting_ = model.phi_;
         }
         public override double shortRate(double t, Vector v)
         {
            return base.shortRate(t, v) + fitting_.value(t);
         }
      }
      #endregion

      #region IAffineModel implémentation
      // pas sur, initialement PhiKshi(now, maturity)
      private double Phiksi(double u, double v)
      {
         double P0u = this.TermStructureInitialForwardRate(u);
         double P0v = this.TermStructureInitialForwardRate(v);
         double discountu = ((CIR2)this).Discount(u);
         double discountv = ((CIR2)this).Discount(v);
         return P0v * discountu / (P0u * discountv);
      }
      public override double A(double t, double T)
      {
         return base.A(t, T) * Phiksi(t,T);
      }
      public override double DiscountBondOption(Option.Type type, double strike, double maturity, double bondMaturity)
      {
         double nominal = 1;
         double ZBC =  nominal * Phiksi(0.0,bondMaturity) * base.DiscountBondOption(type, strike/Phiksi(maturity,bondMaturity), maturity, bondMaturity);
         switch(type)
         {
            case Option.Type.Call:
               return ZBC;
            case Option.Type.Put:
               return ZBC - nominal * ((CIR2)this).Discount(bondMaturity) + strike * ((CIR2)this).Discount(maturity);
               // On utilise le facteur discount du CIR2 qui utilise le A(t,T) du CIR2 et pas dur CIR2PlusPlus
         }
         throw new NotSupportedException("unsupported option type");  
      }
      #endregion

   }
}
