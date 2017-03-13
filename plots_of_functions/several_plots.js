function nondim(array1, adim){
  var output = [];
  for (var i = 0; i < array1.length; i++) {
    output.push(array1[i]/adim);
  }
  return output
};




function plot_cl_alpha_cd(IDDIV) {
  var alpha = createArraySequence(-15,.25,15);
  var cl = [];
  var cd = [];

  for (var i = 0; i < alpha.length; i++) {
   temp = polarAirfoil(alpha[i]);
   cl.push(temp[0]);
   cd.push(temp[1]);
  }


  var trace1 = {
    x: alpha,
    y: cl,
    mode: 'line',
    name: 'Cl'
  };

  var trace2 = {
    x: alpha,
    y: cd,
    mode: 'line',
    name: 'Cd'
  };

  var data = [trace1 , trace2];

  var layout = {
    // title: '$ \\rm{Thrust~coefficient~} C_T \\rm{~as~a~function~of~induction~factor~} a $',
    xaxis: {title: 'Angle of attack',
    dtick: 2,
    range: [-16., 16],
    autorange: false
  },
  yaxis: {
    title: "Cl",
    range: [-1., 1.6],
    dtick: 0.2,
    autorange: false
  }
};

Plotly.newPlot(IDDIV, data, layout);
};




function plot_cl_cd(IDDIV) {
  var alpha = createArraySequence(-15,.25,15);
  var cl = [];
  var cd = [];

  for (var i = 0; i < alpha.length; i++) {
   temp = polarAirfoil(alpha[i]);
   cl.push(temp[0]);
   cd.push(temp[1]);
  }


  var trace1 = {
    x: cd,
    y: cl,
    mode: 'line',
    name: ''
  };

  var data = [trace1];

  var layout = {
    // title: '$ \\rm{Thrust~coefficient~} C_T \\rm{~as~a~function~of~induction~factor~} a $',
    xaxis: {title: 'Cd',
    dtick: 0.002,
    range: [0., 0.02],
    autorange: false
  },
  yaxis: {
    title: "Cl",
    range: [-1., 1.6],
    dtick: 0.2,
    autorange: false
  }
};

Plotly.newPlot(IDDIV, data, layout);
};




function CT_function_a(IDDIV) {

  var a = createArraySequence(0,.01,1);
  var CT = thrust_coefficient_from_induction_Gluert_correction(a);

  var trace1 = {
    x: a,
    y: CT,
    mode: 'line',
    name: 'Scatter'
  };


  var data = [trace1];

  var layout = {
    // title: '$ \\rm{Thrust~coefficient~} C_T \\rm{~as~a~function~of~induction~factor~} a $',
      xaxis: {title: 'Axial induction factor',
      dtick: 0.1,
      range: [0., 1.0],
    autorange: false
  },
    yaxis: {
      title: 'Thrust coefficient',
      range: [0., 2.0],
      dtick: 0.2,
      autorange: false
    }
  };

  Plotly.newPlot(IDDIV, data, layout);
};


function plotBEM(inputx, inputy1, inputy2 , name1, name2, IDDIV) {


  var trace1 = {
    x: inputx,
    y: inputy1,
    mode: 'line',
    name: name1
  };

  var trace2 = {
    x: inputx,
    y: inputy2,
    mode: 'line',
    name: name2
  };

  var data = [trace1, trace2];

  var layout = {
    // title: '$ \\rm{Thrust~coefficient~} C_T \\rm{~as~a~function~of~induction~factor~} a $',
    xaxis: {title: 'r/R',
    dtick: 0.1,
    range: [0., 1.0],
    autorange: false
  },
  yaxis: {
    title: "",
    range: [0., 1.0],
    dtick: 0.2,
    autorange: false
  }
};

Plotly.newPlot(IDDIV, data, layout);
};




function plotFaxial(inputx, inputy1, inputy2 , name1, name2, IDDIV) {


  var trace1 = {
    x: inputx,
    y: inputy1,
    mode: 'line',
    name: name1
  };

  var trace2 = {
    x: inputx,
    y: inputy2,
    mode: 'line',
    name: name2
  };

  var data = [trace1, trace2];

  var layout = {
    // title: '$ \\rm{Thrust~coefficient~} C_T \\rm{~as~a~function~of~induction~factor~} a $',
    xaxis: {title: 'r/R',
    dtick: 0.1,
    range: [0., 1.0],
    autorange: false
  },
  yaxis: {
    title: "",
    range: [0., 2.0],
    dtick: .1,
    autorange: false
  }
};

Plotly.newPlot(IDDIV, data, layout);
};



function plotGamma(inputx, inputy1, name1, IDDIV) {


  var trace1 = {
    x: inputx,
    y: inputy1,
    mode: 'line',
    name: name1
  };

  // var trace2 = {
  //   x: inputx,
  //   y: inputy2,
  //   mode: 'line',
  //   name: name2
  // };

  var data = [trace1];

  var layout = {
    // title: '$ \\rm{Thrust~coefficient~} C_T \\rm{~as~a~function~of~induction~factor~} a $',
    xaxis: {title: 'r/R',
    dtick: 0.1,
    range: [0., 1.0],
    autorange: false
  },
  yaxis: {
    title: "",
    range: [0., 1.2],
    dtick: .1,
    autorange: false
  }
};

Plotly.newPlot(IDDIV, data, layout);
};






function clearBox(elementID)
{
    document.getElementById(elementID).innerHTML = "";
}

function plot_Prandtl_correction(IDDIV){
  var temp;
  var a = createArraySequence(0.01,.00025,1);
  var FX = [];
// console.log(a)

  for (var i = 0; i < a.length; i++) {
    temp=calculatePrandtlTipRootCorrection(a[i], 0.01, 1, 6, 3, 0.1);
    FX.push(temp.Ftotal);
  };
  //  console.log(temp)
// console.log(FX)
  var trace1 = {
    x: a,
    y: FX,
    mode: 'line',
    name: 'Scatter'
  };


  var data = [trace1];

  var layout = {
    // title: '$ \\rm{Thrust~coefficient~} C_T \\rm{~as~a~function~of~induction~factor~} a $',
    xaxis: {title: 'r/R',
    dtick: 0.1,
    range: [0., 1.0],
    autorange: false
  },
  yaxis: {
    title: "Prandtl's correction factor",
    range: [0., 2.0],
    dtick: 0.2,
    autorange: false
  }
};

Plotly.newPlot(IDDIV, data, layout);
};




function plot_Prandtl_correction2(IDDIV){
  var temp;
  var a = createArraySequence(0.01,.00025,1);
  var FX = [];
// console.log(a)

  for (var i = 0; i < a.length; i++) {
    temp=calculatePrandtlTipRootCorrection(a[i], 0.01, 1, 6, 3, 0.1);
    FX.push(temp.Ftotal*1.5);
  }
  //  console.log(temp)
// console.log(FX)
  var trace1 = {
    x: a,
    y: FX,
    mode: 'line',
    name: 'Scatter'
  };


  var data = [trace1];

  var layout = {
    // title: '$ \\rm{Thrust~coefficient~} C_T \\rm{~as~a~function~of~induction~factor~} a $',
    xaxis: {title: 'r/R',
    dtick: 0.1,
    range: [0., 1.0],
    autorange: false
  },
  yaxis: {
    title: "Prandtl's correction factor",
    range: [0., 2.0],
    dtick: 0.2,
    autorange: false
  }
};

Plotly.newPlot(IDDIV, data, layout);
};
