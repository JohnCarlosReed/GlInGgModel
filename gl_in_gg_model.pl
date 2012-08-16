# --------------------------------------------------------------------
# FILENAME: gl_in_gg_model.pl
# 
# AUTHOR:   John Reed
#
# PURPOSE:  This is a perl version mock-up of Professor Cant's 
#           glucose/insulin/glucagon model.  The purpose is for
#           me to build a working model so that I may understand the
#           arithmetic involved. 
#           I chose to do this in Perl first because I can easily
#           build a quick and dirty version of what I need to know.
#           There are still some bugs in this version that need to 
#           be panned out, they are commented.
#
# NOTE:     The fortran version is attached as a comment to the 
#           bottom of this script.
#
# DATE:     06.21.01   
# --------------------------------------------------------------------
use strict;

# ---  GLOBALS  ------------------
my $tstp    = 20;
my $vol     = 8;
my $b0abs   = 1;
my $b1abs   = 0.03;
my $cabs    = 0.015;
my $fedGl   = 0;
my $iGl     = 40.0;
my $infGl   = 0.0;
my $JGgUGl  = 125;
my $JInGNG  = 40.0;
my $KGgGNG  = 125; 
my $KGlUGl  = 11.0;
my $KInUGl  = 40.0;
my $tabson  = 120;
my $VGl     = 7.072;
my $VGNG    = 2.21;
my $fGNGLc  = 0.10;
my $iHotGl  = 0.0000000001;
my $iHotLc  = 0.0000000001;
my $Lc      = 10.0;
my $SAabsGl = 0.0;
my $SAGNG   = 0.0;
my $SAinfGl = 0.0;
my $basalIn = 0.0075;
my $expIn   = 5.0;
my $iIn     = 64.0;
my $infIn   = 0;
my $kIn     = 1.0;
my $KGlPIn  = 8.0;

my $VIn     = 91.8;
my $iInsig  = 8.0;
my $kPInsig = 1.0;
my $kUInsig = 1.0;
my $expGg   = 4.0;
my $iGg     = 640;
my $infGg   = 1000000; #0;
my $JGlPGg  = 4.0;
my $kGg     = 0.625;
my $VGg     = 172.07;
my $iGgsig  = 80;
my $kPGgsig = 0.5;
my $kUGgsig = 0.5;
my ( $GNG, $Insig, $Ggsig, $UGl );
my ( $PGg, $UGg, $cGl, $cGg );
my ( $PIn, $UIn, $cIn );
my ( $PInsig, $UInsig, $PGgsig,$UGgsig );
my ( $absGl, $SAGl );

my ( $dGldt1,    $dGldt2,    $dGldt3,    $dGldt4 );
my ( $dHotGldt1, $dHotGldt2, $dHotGldt3, $dHotGldt4 );
my ( $dIndt1,    $dIndt2,    $dIndt3,    $dIndt4 ); 
my ( $dInsigdt1, $dInsigdt2, $dInsigdt3, $dInsigdt4 );
my ( $dGgdt1,    $dGgdt2,    $dGgdt3,    $dGgdt4 ); 
my ( $dGgsigdt1, $dGgsigdt2, $dGgsigdt3, $dGgsigdt4 );

my ( $Gl, $HotGl, $In, $Gg );  # I might be able to get rid of these

my $dt        = 0.005;
my $days      = 7;             # I made these up, need to ask Cant 
my $intervals = $days/$dt;     # because I dont speak Fortran.  


# ---------------------------------------------------------------------------
# Set the Initial states
#
$Ggsig = $iGgsig;
$Insig = $iInsig;
$cGl   = $iGl/$vol;            
$cIn   = $iIn/$vol;
$cGg   = $iGg/$vol;
$SAGl  = $iHotGl/$iGl;         


# -------------------------------------------------------------------------------
# Foreach time interval, do the four Runge-Kutta calculations 
#
for ( 0..$intervals ){


  # -- K1 ------------------------------------------------------------------------   
  evaluate_fluxes();

  $dGldt1    = $absGl + $infGl + $GNG - $UGl;
  $dHotGldt1 = $absGl*$SAabsGl + $infGl*$SAinfGl + $GNG*$SAGNG - $UGl*$SAGl;
  $dIndt1    = $infIn + $PIn - $UIn;
  $dInsigdt1 = $PInsig - $UInsig;
  $dGgdt1    = $infGg + $PGg - $UGg;
  $dGgsigdt1 = $PGgsig - $UGgsig;
 
  my $kGl1    = $iGl    + $dt * ( $dGldt1    / 2 );
  my $kHotGl1 = $iHotGl + $dt * ( $dHotGldt1 / 2 );  # THIS LOOKS WRONG
  my $kIn1    = $iIn    + $dt * ( $dIndt1    / 2 );
  my $kInsig1 = $iInsig + $dt * ( $dInsigdt1 / 2 );
  my $kGg1    = $iGg    + $dt * ( $dGgdt1    / 2 );
  my $kGgsig1 = $iGgsig + $dt * ( $dGgsigdt1 / 2 );


  $Gl    = $kGl1;
  $cGl   = $Gl/$vol;

  $HotGl = $kHotGl1;
  $SAGl  = $HotGl/$Gl;

  $In    = $kIn1;    
  $cIn   = $In/$vol;

  $Insig = $kInsig1;

  $Gg    = $kGg1;
  $cGg   = $Gg/$vol;  
   
  $Ggsig = $kGgsig1; 

  # ---------------------------------------------------------------------------------


  # -- K2 ---------------------------------------------------------------------------

  evaluate_fluxes();

  $dGldt2    = $absGl + $infGl + $GNG - $UGl;
  $dHotGldt2 = $absGl*$SAabsGl + $infGl*$SAinfGl + $GNG*$SAGNG - $UGl*$SAGl;
  $dIndt2    = $infIn + $PIn - $UIn;
  $dInsigdt2 = $PInsig - $UInsig;
  $dGgdt2    = $infGg + $PGg - $UGg;
  $dGgsigdt2 = $PGgsig - $UGgsig;

  my $kGl2    = $iGl    + $dt * ( $dGldt2    / 2 );
  my $kHotGl2 = $iHotGl + $dt * ( $dHotGldt2 / 2 );
  my $kIn2    = $iIn    + $dt * ( $dIndt2    / 2 );

  my $kInsig2 = $iInsig + $dt * ( $dInsigdt2 / 2 );
  my $kGg2    = $iGg    + $dt * ( $dGgdt2    / 2 );
  my $kGgsig2 = $iGgsig + $dt * ( $dGgsigdt2 / 2 );


  $Gl  = $kGl2;
  $cGl = $Gl/$vol;

  $HotGl = $kHotGl2;
  $SAGl  = $HotGl/$Gl;

  $In  = $kIn2;
  $cIn = $In/$vol;

  $Insig = $kInsig2;

  $Gg  = $kGg2;
  $cGg = $Gg/$vol;

  $Ggsig = $kGgsig2; 
  # ---------------------------------------------------------------------------------
  

  # -- K3 ---------------------------------------------------------------------------
  evaluate_fluxes();

  $dGldt3    = $absGl + $infGl + $GNG - $UGl; 

  $dHotGldt3 = $absGl*$SAabsGl + $infGl*$SAinfGl + $GNG*$SAGNG - $UGl*$SAGl;
  $dIndt3    = $infIn + $PIn - $UIn;
  $dInsigdt3 = $PInsig - $UInsig;
  $dGgdt3    = $infGg + $PGg - $UGg;
  $dGgsigdt3 = $PGgsig - $UGgsig;

  my $kGl3    = $iGl    + $dt * ( $dGldt3 );
  my $kHotGl3 = $iHotGl + $dt * ( $dHotGldt3 );
  my $kIn3    = $iIn    + $dt * ( $dIndt3 );
  my $kInsig3 = $iInsig + $dt * ( $dInsigdt3 );
  my $kGg3    = $iGg    + $dt * ( $dGgdt3 );
  my $kGgsig3 = $iGgsig + $dt * ( $dGgsigdt3 );



  $Gl  = $kGl3;
  $cGl = $Gl/$vol;

  $HotGl = $kHotGl3;
  $SAGl  = $HotGl/$Gl;

  $In  = $kIn3;
  $cIn = $In/$vol;

  $Insig = $kInsig3;

  $Gg  = $kGg3;
  $cGg = $Gg/$vol;

  $Ggsig = $kGgsig3; 
  # ---------------------------------------------------------------------------------


  # -- K4 ---------------------------------------------------------------------------
  evaluate_fluxes();

  $dGldt4    = $absGl + $infGl + $GNG - $UGl;
  $dHotGldt4 = $absGl*$SAabsGl + $infGl*$SAinfGl + $GNG*$SAGNG - $UGl*$SAGl;
  $dIndt4    = $infIn + $PIn - $UIn;
  $dInsigdt4 = $PInsig - $UInsig;
  $dGgdt4    = $infGg + $PGg - $UGg;
  $dGgsigdt4 = $PGgsig - $UGgsig;

  my $kGl4    = $iGl    + ($dt/6) * ( $dGldt1    + 2*$dGldt2    + 2*$dGldt3    + $dGldt4    );
  my $kHotGl4 = $iHotGl + ($dt/6) * ( $dHotGldt1 + 2*$dHotGldt2 + 2*$dHotGldt3 + $dHotGldt4 );
  my $kIn4    = $iIn    + ($dt/6) * ( $dIndt1    + 2*$dIndt2    + 2*$dIndt3    + $dIndt4    );
  my $kInsig4 = $iInsig + ($dt/6) * ( $dInsigdt1 + 2*$dInsigdt2 + 2*$dInsigdt3 + $dInsigdt4 );
  my $kGg4    = $iGg    + ($dt/6) * ( $dGgdt1    + 2*$dGgdt2    + 2*$dGgdt3    + $dGgdt4    );
  my $kGgsig4 = $iGgsig + ($dt/6) * ( $dGgsigdt1 + 2*$dGgsigdt2 + 2*$dGgsigdt3 + $dGgsigdt4 );


  $Gl  = $kGl4;
  $cGl = $Gl/$vol;

  $HotGl = $kHotGl4;
  $SAGl  = $HotGl/$Gl;

  $In  = $kIn4;
  $cIn = $In/$vol;

  $Insig = $kInsig4;

  $Gg  = $kGg4;
  $cGg = $Gg/$vol;

  $Ggsig = $kGgsig4; 
  # ---------------------------------------------------------------------------------

  # reset initials
  $iGl = $Gl;
  $iHotGl = $HotGl;
  # SAGl here? 
  $iIn    = $In;
  $iInsig = $Insig;
  $iGg    = $Gg;
  $iGgsig = $Ggsig;


#printf( "Gl:%.7f HotGl:%.7f SAGl:%.7f In:%.7f Insig:%.2f Gg:%.2f \n", 
#          $Gl,    $HotGl,    $SAGl,    $In,    $Insig,    $Gg );


  printf( "cGl:%.7f cIn:%.7f cGg:%.7f SAGl:%.7f GNG:%.2f UGl:%.2f \n",
          $cGl,    $cIn,    $cGg,    $SAGl,    $GNG,    $UGl );


} # ends interval loop



# -----------------------------------------------------------------------------------
# This is the core of the model. All of the above code is calculus required to 
# determine the change in these rates
# 
sub evaluate_fluxes{

  $absGl = 0;  # need to put that fortran PROCEDURAL here.

  $GNG    = $VGNG/(1+$KGgGNG/$Ggsig+$Insig/$JInGNG);
  $UGl    = $VGl/(1+$KGlUGl/$cGl+$KInUGl/$Insig+$Ggsig/$JGgUGl);
  $PIn    = $VIn/(1+($KGlPIn/$cGl)**$expIn) + $basalIn;   
  $UIn    = $kIn*$cIn;
  $PInsig = $kPInsig*$cIn;
  $UInsig = $kUInsig*$Insig;
  $PGg    = $VGg/(1+($cGl/$JGlPGg)**$expGg);              
  $UGg    = $kGg*$cGg;
  $PGgsig = $kPGgsig*$cGg;
  $UGgsig = $kUGgsig*$Ggsig;
}

# -- END --------------------------------------------------------------------------------




# ----------------------------------------------------------------------------------------
# FORTRAN VERSION
# ----------------------------------------------------------------------------------------
#! glucose/insulin/glucagon dynamic model
#! with In and Gg signal pools
#! labelled glc
#PROGRAM minmod4

#INITIAL
#       ! units are in dpm, L, min, mmol, mU (for In), ng (for Gg)
#       CONSTANT tstp=20.0, vol=8.0
#       LOGICAL abson
#END $ 'of initial'
#DYNAMIC
#       ALGORITHM IALG=5

#       NSTEPS NSTP=1
#       CONSTANT CINT=1.0
#       MAXTERVAL MAXT=0.005
#    DERIVATIVE
#       TERMT (T.GE.TSTP)
#
#!**** GLUCOSE *****
#CONSTANT abson=.false., b0abs=1.0, b1abs=0.03, cabs=0.015, fedGl=0
#CONSTANT iGl=40.0, infGl=0.0, JGgUGl=125, JInGNG=40.0
#CONSTANT KGgGNG=125, KGlUGl=11.0, KInUGl=40.0, tabson=120
#CONSTANT VGl=7.072, VGNG=2.21

#       dGldt = absGl + infGl + GNG - UGl
#       PROCEDURAL(absGl=abson)
#               IF (abson.eq..true.) THEN
#                       Aabs = fedGl/(1/cabs-1/(b1abs+cabs))
#                       absGl = abson*(Aabs*(b0abs-exp(-b1abs*(t-tabson)))*...
#                               exp(-cabs*(t-tabson)))
#                       ! absGl in mmol/min, fedGl in mmol (per meal/drink)
#               ELSE
#                       absGl=0.0
#               ENDIF
#       END ! of procedural
#       GNG = VGNG/(1+KGgGNG/Ggsig+Insig/JInGNG)
#       UGl = VGl/(1+KGlUGl/cGl+KInUGl/Insig+Ggsig/JGgUGl)
#       Gl = INTEG(dGldt, iGl)
#       cGl = Gl/vol
#
#!**** LABELLED GLUCOSE ****
#CONSTANT fGNGLc=0.10, iHotGl=1.0e-10, iHotLc=1.0e-10, Lc = 10.0
#CONSTANT SAabsGl=0.0, SAGNG=0.0, SAinfGl=0.0
#       dHotGldt = absGl*SAabsGl + infGl*SAinfGl + GNG*SAGNG -  UGl*SAGl
#       HotGl = INTEG(dHotGldt, iHotGl)
#       SAGl = HotGl/Gl
#
#!**** INSULIN ****
#CONSTANT basalIn=0.0075, expIn=5.0, iIn=64.0, infIn=0, kIn=1.0
#CONSTANT KGlPIn=8.0, VIn=91.8
#       dIndt = infIn + PIn - UIn
#       PIn = VIn/(1+(KGlPIn/cGl)**expIn) + basalIn
#       UIn = kIn*cIn
#       In = INTEG(dIndt, iIn)
#       cIn = In/vol
#
#CONSTANT iInsig=8.0, kPInsig=1.0, kUInsig=1.0
#       dInsigdt = PInsig - UInsig
#       PInsig = kPInsig*cIn
#       UInsig = kUInsig*Insig
#       Insig = INTEG(dInsigdt, iInsig)
#
#!**** GLUCAGON ****
#CONSTANT expGg=4.0, iGg=640, infGg=0, JGlPGg=4.0, kGg=0.625, VGg=172.07
#       dGgdt = infGg + PGg - UGg
#       PGg = VGg/(1+(cGl/JGlPGg)**expGg)
#       UGg = kGg*cGg
#       Gg = INTEG(dGgdt, iGg)
#       cGg = Gg/vol
#
#CONSTANT iGgsig=80, kPGgsig=0.5, kUGgsig=0.5
#       dGgsigdt = PGgsig - UGgsig
#       PGgsig = kPGgsig*cGg
#       UGgsig = kUGgsig*Ggsig
#       Ggsig = INTEG(dGgsigdt, iGgsig)
#
#    END $ 'of derivative'
#END $ 'of dynamic'
#END $ 'of program'
#