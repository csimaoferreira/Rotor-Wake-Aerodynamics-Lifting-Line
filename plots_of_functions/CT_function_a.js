// plot of CT as a function of a

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

Plotly.newPlot('chartContainer', data, layout);
