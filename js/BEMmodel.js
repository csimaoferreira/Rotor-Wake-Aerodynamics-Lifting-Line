// BEM fuctions


function thrust_coefficient_from_induction(a) {
  // calculates thrust coefficient as a function of induction factor a
  var CT=[];

  for (var i = 0; i < a.length; i++) {
    temp = 4*a[i]*(1-a[i]);
    CT.push(temp);
  };
  return CT;
};



function thrust_coefficient_from_induction_Gluert_correction(a) {
  // calculates thrust coefficient as a function of induction factor a
  var CT=[];
  var CT1=1.816;
  var a1=1-Math.sqrt(CT1)/2;

  for (var i = 0; i < a.length; i++) {

    if (a[i]<a1) {
      temp = 4*a[i]*(1-a[i]);
    } else {
      temp = CT1-4*(Math.sqrt(CT1)-1)*(1-a[i]);
    };

    CT.push(temp);
  };
  return CT;
};


function induction_from_thrust_coefficient_Gluert_correction(CT) {
  // calculates induction factor a as a function of thrust coefficient
  var a=[];
  var CT1=1.816;
  var CT2=2*Math.sqrt(CT1)-CT1;

  for (var i = 0; i < CT.length; i++) {

    if (CT[i]<CT2) {
      temp = 0.5-0.5*Math.sqrt(1-CT[i]);
    } else {
      temp = 1 + (CT[i]-CT1)/(4*(Math.sqrt(CT1)-1));
    };

    a.push(temp);
  };
  return a;
};


function createArraySequence(valuestart,deltavalue,valueend){
  var data = [];
  var temp = valuestart;
  while (temp <= valueend) {
    data.push(temp);
    temp =temp+deltavalue;
  }
  return data;                // Function returns data
};

function interpolateArray(xarray,yarray,xnew){
  // interpolate 1D array
  var i1;
  var i2;
  var result;
  for (var i = 0; i < (xarray.length-1); i++) {
    if (xarray[i]<=xnew) {
      if (xarray[i+1]>=xnew) {
        i1=i;
        i2=i+1;
      };
    };
  };
  result=(xnew-xarray[i1])/(xarray[i2]-xarray[i1])*(yarray[i2]-yarray[i1])+yarray[i1];
  return result;
};




function geoBlade(r_R) {
  var radial = [0, 0.3, .5, .8, 1];
  // var chorddist = [.05, .04, .03, .02, .015];
  // var twistdist = [-12, -8, -5, -4, 0.];
  var pitch = 2;
  var chord = 3*(1-r_R)+1;
  var twist = -14*(1-r_R);
  var result =[chord , twist + pitch];
  return result;
};


function loadBladeElement(Vnorm, Vtan, r_R){
  var Vmag2 = (Math.pow(Vnorm,2) + Math.pow(Vtan,2)  );
  var InflowAngle = Math.atan(Vnorm/Vtan);
  // get chord and twist
  // console.log('inflow angle ' + InflowAngle)
  var temp= geoBlade(r_R);
  var chord = temp[0];
  var twist = temp[1];
  // console.log('twist ' + twist)
  var alpha = twist + InflowAngle*180/Math.PI;
  // console.log('alpha ' + alpha)
  temp = polarAirfoil(alpha);
  var cl = temp[0];
  var cd = 0*temp[1];
  // console.log('cl and cd ' + cl +' ' +cd)
  var Lift = 0.5*Vmag2*cl*chord;
  var Drag = 0.5*Vmag2*cd*chord;
  var Fnorm = Lift*Math.cos(InflowAngle)+Drag*Math.sin(InflowAngle);
  var Ftan = Lift*Math.sin(InflowAngle)-Drag*Math.cos(InflowAngle);
  var Gamma = 0.5*Math.sqrt(Vmag2)*cl*chord ; //
  var result = [Fnorm , Ftan, Gamma];
  // console.log('Result ' + result)

  return result;
};




function calculateCT_CProtor_CPflow(a,aline,Fnorm,Ftan, Uinf, r_Rarray, Omega, Radius, NBlades){
  // calculates the performance of the rotor, retunring CT and CPs
  var CTrotor=0; // trhust coefficient
  var r_R_temp=0; // radial position for evaluation
  var drtemp=0; // delta r
  var CProtor = 0;
  var CPflow = 0;
  var temp
  var Utan =0
  var Unorm =0;
  for (var i = 0; i < (r_Rarray.length-1); i++) {
    r_R_temp = (r_Rarray[i]+r_Rarray[i+1])/2;
    Utan = r_R_temp*Radius*Omega*aline[i];
    Unorm = Uinf*(1-a[i]);
    drtemp = (-r_Rarray[i]+r_Rarray[i+1]);
    CTrotor+=(drtemp*Fnorm[i]*NBlades)/(0.5*Uinf*Uinf*Math.PI*Radius);
    CProtor+=(drtemp*Ftan[i]*r_R_temp*Omega*NBlades)/(0.5*Uinf*Uinf*Uinf*Math.PI);
    // CPflow+=(drtemp*(Fnorm[i]*Unorm - Ftan[i]*Utan )*NBlades)/(0.5*Uinf*Uinf*Uinf*Math.PI*Radius);



  }
  results ={CTrotor: CTrotor, CProtor: CProtor, CPflow: CPflow };
  return results;
};








function solveBEMmodel(Uinf, r_Rarray, Omega, Radius, NBlades){
  // calculates the solution of BEM over the entire rotor, returning
  // arrays of axial induction factor a, azimuthal induction factor a'
  var a_temp=[]; //axial induction factor a
  var aline_temp=[]; //azimuthal induction factor a'
  var r_R_temp=[]; // radial position for evaluation
  var Fnorm = [];
  var Ftan = [];
  var temp;
  var Gamma = []; // circulation
  for (var i = 0; i < (r_Rarray.length-1); i++) {
    temp=solveStreamtube(Uinf, r_Rarray[i], r_Rarray[i+1],r_Rarray[0], r_Rarray[r_Rarray.length-1] , Omega, Radius, NBlades )
    a_temp.push(temp[0]);
    aline_temp.push(temp[1]);
    r_R_temp.push(temp[2]);
    Fnorm.push(temp[3]);
    Ftan.push(temp[4]);
    Gamma.push(temp[5]);


  }
  results ={a: a_temp , aline: aline_temp, r_R: r_R_temp, Fnorm: Fnorm, Ftan: Ftan , Gamma: Gamma};
  return results;
};


function solveStreamtube(Uinf, r1_R, r2_R, rootradius_R, tipradius_R , Omega, Radius, NBlades ){
  // solve balance of momentum between blade element load and loading in the streamtube
  // input variables:
  //     Uinf - wind speed at infinity
  //     r1_R,r2_R - edges of blade element, in fraction of Radius ;
  //     rootradius_R, tipradius_R - location of blade root and tip, in fraction of Radius ;
  //     Radius is the rotor radius
  //     Omega -rotational velocity
  //     NBlades - number of blades in rotor

  // initialize properties of the blade element, variables for output and induction factors
  var r_R = (r1_R+r2_R)/2; //centroide
  var Area = Math.PI*(Math.pow(r2_R*Radius,2)-Math.pow(r1_R*Radius,2)); //  area streamtube
  var a = 0.3; // axial induction factor
  var anew; // temp new axial induction factor
  var aline = 0.; // tangential induction factor
  var Urotor; // axial velocity at rotor
  var Utan; // tangential velocity at rotor
  var loads; // normal and tangential loads 2D
  var load3D = [0 , 0]; // normal and tangential loads 3D
  var CT; //thrust coefficient at streamtube
  var Prandtl; // Prandtl tip correction

  // iteration cycle
  var Niterations =100; // maximum number of iterations
  var Erroriterations =0.00001; // error limit for iteration rpocess, in absolute value of induction
  for (var i = 0; i < Niterations; i++) {

    ///////////////////////////////////////////////////////////////////////
    // this is the block "Calculate velocity and loads at blade element"
    ///////////////////////////////////////////////////////////////////////
    Urotor = Uinf*(1-a); // axial velocity at rotor
    Utan = (1+aline)*Omega*r_R*Radius; // tangential velocity at rotor
    // calculate loads in blade segment in 2D (N/m)
    loads = loadBladeElement(Urotor, Utan, r_R);
    load3D[0] =loads[0]*Radius*(r2_R-r1_R)*NBlades; // 3D force in axial direction
    load3D[1] =loads[1]*Radius*(r2_R-r1_R)*NBlades; // 3D force in azimuthal/tangential direction (not used here)
    ///////////////////////////////////////////////////////////////////////
    //the block "Calculate velocity and loads at blade element" is done
    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////
    // this is the block "Calculate new estimate of axial and azimuthal induction"
    ///////////////////////////////////////////////////////////////////////
    // calculate thrust coefficient at the streamtube
    CT = load3D[0]/(0.5*Area*Math.pow(Uinf,2));
    // calculate new axial induction, accounting for Glauert's correction
    anew = induction_from_thrust_coefficient_Gluert_correction([CT]);
    // correct new axial induction with Prandtl's correction
    Prandtl=calculatePrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, Omega*Radius/Uinf, NBlades, anew);
    if (Prandtl.Ftotal < 0.0001) { Prandtl.Ftotal = 0.0001; } // avoid divide by zero
    anew = anew/Prandtl.Ftotal; // correct estimate of axial induction
    a = 0.75*a+0.25*anew; // for improving convergence, weigh current and previous iteration of axial induction
    // calculate aximuthal induction
    aline = loads[1]*NBlades/(2*Math.PI*Uinf*(1-a)*Omega*2*Math.pow(r_R*Radius,2));
    aline =aline/Prandtl.Ftotal; // correct estimate of azimuthal induction with Prandtl's correction
    ///////////////////////////////////////////////////////////////////////////
    // end of the block "Calculate new estimate of axial and azimuthal induction"
    ///////////////////////////////////////////////////////////////////////

    // test convergence of solution, by checking convergence of axial induction
    if (Math.abs(a-anew) < Erroriterations) {
      // console.log(i);
      i=Niterations; // converged solution, this is the last iteration
    }

  };

  // we have reached a solution or the maximum number of iterations
  // returns axial induction factor a, azimuthal induction factor a',
  // and radial position of evaluations and loads
  return [a , aline, r_R, loads[0] , loads[1], loads[2]];
};


function calculatePrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, TSR, NBlades, axial_induction) {
// applies Prandtl's tip and root correction to the induction vector aind
var mu; // temp variabl
var Ftip; // tip correction
var Froot; // tip correction
var result;
var temp1;
    temp1 = -NBlades/2*(tipradius_R-r_R)/r_R*Math.sqrt( 1+ (Math.pow(TSR*r_R,2))/(Math.pow(1-axial_induction,2))   )
    Ftip = 2/Math.PI*Math.acos(Math.exp(temp1));
    if (isNaN(Ftip)) {
      Ftip = 0;
    };
    temp1 = NBlades/2*(rootradius_R-r_R)/r_R*Math.sqrt( 1+ (Math.pow(TSR*r_R,2))/(Math.pow(1-axial_induction,2))   )
    Froot = 2/Math.PI*Math.acos(Math.exp(temp1));
    if (isNaN(Froot)) {
      Froot = 0;
    };
    result= {Ftotal: (Froot*Ftip), Ftip: Ftip ,Froot: Froot};
    return result;
}




function polarAirfoil(alpha) {
  var polar_alpha = [-180, -16.062 , 	-15.506 , 	-15.064 , 	-14.589 , 	-14.109 , 	-13.698 ,
    -13.237 , 	-12.745 , 	-12.268 , 	-11.748 , 	-11.183 , 	-10.768 , 	-10.231 , 	-9.743 ,
    -9.223 , 	-8.209 , 	-7.187 , 	-6.162 , 	-5.143 , 	-4.127 , 	-3.106 , 	-2.073 , 	-1.04 ,
    .017 , 	1.025 , 	2.042 , 	3.096 , 	4.114 , 	5.126 , 	6.163 , 	7.189 , 	7.713 ,
    8.216 , 	8.734 , 	9.251 , 	9.558 , 	9.771 , 	10.269 , 	10.757 , 	11.257 , 	11.761 ,
    12.239 , 	13.224 , 	14.234 , 	15.227 , 	16.208 , 	17.201 , 	18.2 , 	19.201 , 	20.189 ,
    21.179 , 	22.162 , 	23.144 , 	23.547 , 	24.062 , 	25.056 , 	26.06 , 	27.059 ,
    28.062 , 	29.062 , 	30.056 , 180];
    var polar_cl = [-0, -.425 , 	-.43 , 	-.461 , 	-.496 , 	-.567 , 	-.682 , 	-.719 , 	-.755 ,
      -.77 , 	-.774 , 	-.77 , 	-.756 , 	-.736 , 	-.711 , 	-.68 , 	-.612 , 	-.529 , 	-.434 ,
      -.337 , 	-.237 , 	-.127 , 	-.014 , 	.102 , 	.217 , 	.33 , 	.445 , 	.571 , 	.683 ,
      .792 , 	.902 , 	1.014 , 	1.068 , 	1.12 , 	1.168 , 	1.207 , 	1.227 , 	1.229 ,
      1.215 , 	1.192 , 	1.173 , 	1.148 , 	1.126 , 	1.103 , 	1.093 , 	1.077 ,
      1.06 , 	1.038 , 	1.047 , 	1.043 , 	1.04 , 	1.029 , 	1.023 , 	1.007 , 	.829 ,
      .854 , 	.877 , 	.89 , 	.953 , 	.976 , 	1.016 , 	1.041   , 0];
      var polar_cd = [0.1, .225 , 	.217 , 	.211 , 	.208 , 	.201 , 	.077 , 	.05 , 	.039 ,
        .033 , 	.029 , 	.025 , 	.023 , 	.02 , 	.018 , 	.017 , 	.014 , 	.013 , 	.011 , 	.01 , 	.008 ,
        .008 , 	.008 , 	.007 , 	.007 , 	.008 , 	.008 , 	.008 , 	.008 , 	.009 , 	.009 , 	.009 , 	.009 ,
        .01 , 	.01 , 	.01 , 	.012 , 	.013 , 	.019 , 	.025 , 	.036 , 	.045 , 	.053 , 	.066 , 	.076 , 	.089 , 	.104 ,
        .117 , 	.132 , 	.152 , 	.172 , 	.198 , 	.228 , 	.261 , 	.416 , 	.436 , 	.466 , 	.491 , 	.541 , 	.574 , 	.616 , 	.652
        , 0.1];
        var cl = interpolateArray(polar_alpha,polar_cl,alpha);
        var cd = interpolateArray(polar_alpha,polar_cd,alpha);
        var result =[cl , cd];
        return result;
      };
