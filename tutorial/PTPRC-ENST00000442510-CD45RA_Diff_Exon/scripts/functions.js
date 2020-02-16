function fileDownload(blob,name){
  if(window.navigator.msSaveBlob){
    window.navigator.msSaveBlob(blob, name);
  }else{
    var reader = new FileReader();
    reader.onload = function (event) {
      var save = document.createElement('a');
      save.href = event.target.result;
      save.target = '_blank';
      save.download = name;
      var clicEvent = new MouseEvent('click', {
        'view': window,
        'bubbles': true,
        'cancelable': true
      });
      save.dispatchEvent(clicEvent);
      (window.URL || window.webkitURL).revokeObjectURL(save.href);
    };
    reader.readAsDataURL(blob);
  }
}

function svgDownload(){
  var svg = d3.select("svg");
  svg.select(".buttons").style("display", "none");
  var svgString = new XMLSerializer().serializeToString(svg.node());
  svg.select(".buttons").style("display", null);
  var blob = new Blob([svgString], {type: 'image/svg+xml;charset=utf-8'});
  fileDownload(blob, d3.select("head>title").text()+'.svg');
}

function getCoords(elem, svgDoc) {

  var offset = svgDoc.getBoundingClientRect();

  var bbox = elem.getBBox(),
    x = bbox.x,
    y = bbox.y,
    x2 = x+bbox.width,
    y2 = y+bbox.height;

  var matrix = elem.getScreenCTM();

  return {
    x: (matrix.a * x) + (matrix.c * y) + matrix.e - offset.left,
    y: (matrix.b * x) + (matrix.d * y) + matrix.f - offset.top,
    x2: (matrix.a * x2) + (matrix.c * y2) + matrix.e - offset.left,
    y2: (matrix.b * x2) + (matrix.d * y2) + matrix.f - offset.top
  };
}

function applyOpacity(rgb,alpha,old){
  var blending = function(newC,old){
    return alpha * newC + (1 - alpha) * old;
  }
  return {r: blending(rgb.r,old.r), g: blending(rgb.g,old.g), b: blending(rgb.b,old.b)};
}

function viewport(){
  var e = window,
      a = 'inner';
  if ( !( 'innerWidth' in window ) ){
    a = 'client';
    e = document.documentElement || document.body;
  }
  return { width : e[a+'Width'] , height : e[a+'Height'] }
}

function formatter(d){
  if(typeof d == 'number'){
    var dabs = Math.abs(d);
    if((dabs>0 && dabs<1e-3) || dabs>1e+5)
      d = d.toExponential(3);
    else
      d = d3.round(d,3);
  }
  return d;
}

function displayButtons(){
    var left = typeof margin != 'undefined'?margin.left+width+margin.right:width;
    body.append("button")
        .text("SVG")
        .style({"position":"absolute","top":"35px","left":left+"px"})
        .on('click', svgDownload);
    body.append("button")
        .text("PDF")
        .style({"position":"absolute","top":"60px","left":left+"px"})
        .on('click', svg2pdf);
}

function showWindow(txt){
  var docSize = viewport(),
      bg = body.append("div")
        .style({"position":"fixed","top":0,"left":0,"width":docSize.width+"px","height":docSize.height+"px","background":"rgba(0,0,0,0.8)","z-index":10})

  var win = bg.append("div")
    .attr("class","window")
    .style({"background-color":"#fff","margin":"0 auto","margin-top":(docSize.height/5)+"px","padding":"30px","width":(docSize.width/2)+"px","height":(docSize.height/2)+"px","overflow-y":"auto"})

  if(txt)
    win.text(txt);

  bg.append("div")
    .style({"position":"absolute","top":(12+docSize.height/5)+"px","right":(docSize.width/4)+"px","background":"#fff","width":"17px","height":"17px","border-radius":"13px","text-align":"center","cursor":"pointer"})
    .html("&#x2716;")
    .on("click", function(){ bg.remove() });

  return win;
}

function bioinfoLogo(){
  d3.select("body")
    .append("a")
    .style({"display":"block", "position":"fixed", "right":"10px", "bottom":"10px"})
      .attr("href","http://bioinfo.usal.es")
      .attr("title","Please, cite us: \nBarrios, D. & Prieto, C. D3GB: An Interactive Genome Browser for R, Python, and WordPress. J. Comput. Biol. (2017). doi:10.1089/cmb.2016.0213")
      .html(bioinfoLogoSVG)
}
