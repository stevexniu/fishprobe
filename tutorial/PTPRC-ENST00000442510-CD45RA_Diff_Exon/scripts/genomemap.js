if(typeof genomebrowser == 'undefined'){

  var margin = 0, width = 0, height = 0;

  window.onload = function(){ 
    genomemap();
  }
}

function genomemap(){

    margin = {top: 20, right: 50, bottom: 80, left: 50},
    width = 1000 - margin.left - margin.right,
    height = 580 - margin.top - margin.bottom;

var x = d3.scale.linear()
    .range([10, 26]);

var y = d3.scale.linear()
    .range([0, height]);

d3.select("body").append("div")
    .style({"width":(width + margin.left + margin.right)+"px", "margin": "0 auto"})
    .append("div")
    .style({"padding-left": margin.left+"px", "padding-right": margin.right+"px", "margin-top": "20px"})

if(typeof displaySelector != 'undefined'){
  d3.select("body>div>div")
	.append("a")
	.append("button")
        .attr("disabled", true)
        .text("Go to");
  imgInfo(d3.select("body>div>div"), function(sel){
    sel.append('div').style("padding","10% 20%").html(
           '<p>Click or select a region to view in the genome browser.</p>'+
           infoImages.genomeMap+
           '<p>Then click on <button>Go to</button></p>')
  })
}

var svg = d3.select("body>div").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom);
  svg.append("style")
    .text("text { font: 10px sans-serif;} "+
".mapTrack { fill: #090; }"+
".bands>* {stroke: #000; fill: #FFF;} "+
".acen {fill: #DDD;} "+
".gpos25 {fill: rgb(25%, 25%, 25%);} "+
".gpos33 {fill: rgb(33%, 33%, 33%);} "+
".gpos50 {fill: rgb(50%, 50%, 50%);} "+
".gpos66 {fill: rgb(66%, 66%, 66%);} "+
".gpos75 {fill: rgb(75%, 75%, 75%);} "+
".gpos100 {fill: #000;} "+
".gneg {fill: #FFF;} "+
".gvar {fill: #CCC;} "+
".stalk {fill: #666; stroke: none}");

  svg = svg.append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  var data = JSON.parse(d3.select("#data").text());

  if(typeof chromosomes == 'undefined')
    var chromosomes = JSON.parse(d3.select("#chromosomes").text());

  y.domain([0,data.max]);

  var xMin = 0,
      xMax = 0,
      count = data.chromosomes.length;

  if(data.dataDomain){
    if(data.dataDomain[0]<xMin)
      xMin = data.dataDomain[0];
    if(data.dataDomain[1]>xMax)
      xMax = data.dataDomain[1];
  }

  x.domain([xMin,xMax]).nice();

  if (count>24){
    width = count*(width/24);
    d3.select("svg").attr("width", width + margin.left + margin.right)
  }

data.chromosomes.forEach(function(chr, i){

  svg.append("text")
     .attr("transform", "translate("+(i*width/count)+","+(height+16)+") rotate(45)")
     .text(chr);

    var json = chromosomes[chr.toLowerCase()];

    var gChr = svg.append("g")
      .attr("class", "chromosome")
      .attr("id", chr)
      .attr("transform", "translate("+(i*width/count)+","+height+") rotate(-90)")

    var maxDom = d3.max(json,function(d){ return d[1]; }),
        maxRan = y(maxDom),
        fitScale = d3.scale.linear()
              .range([0, maxRan])
              .domain([0, maxDom]);

    displayChr(gChr, fitScale, json);
    if(typeof displaySelector != 'undefined')       
      displaySelector(gChr, fitScale);

    gChr.selectAll(".mapTrack")
       .data(data.data[chr]?data.data[chr]:[])
     .enter().append("path")
      .attr("class", "mapTrack")
      .attr("d", function(d){
        var x1 = x(0),
            x2 = x(d[2]),
            y1 = y(d[0]),
            y2 = y(d[1]);
        if(y2-y1<1){
          y2 = ((y1+y2)/2)+.5;
          y1 = ((y1+y2)/2)-.5;
        }
        return "M"+y1+","+x1+"L"+y2+","+x1+"L"+y2+","+x2+"L"+y1+","+x2+"Z";
      })
});

displayButtons();
bioinfoLogo();
}

function displayButtons(){
    d3.select("body>div>div").append("button")
	.text("SVG")
	.style("float","right")
	.on("click", svgDownload);
    d3.select("body>div>div").append("button")
	.text("PDF")
	.style("float","right")
	.on("click", svg2pdf);
}

function displayChr(gChr, chrScale, json){

  var gBands = gChr.select("g.bands")

  if(gBands.empty())
    gBands = gChr.append("g")
      .attr("class","bands")

  var bands = gBands.selectAll("path")
        .data(json);
  bands.enter().append("path").append("title");
  bands.exit().remove();
  bands.attr("id", function(d){return (d[2] && d[2]!="")?gChr.attr("id")+"_"+d[2]:null;})
        .attr("class", function(d){return (d[3] && d[3]!="")?d[3]:null;})
        .attr("d", function(d){
            var dPath = "",
		start = chrScale(d[0]),
		end = chrScale(d[1]);	    
            if(d[3] && d[3]=="acen"){
              if(d[2].indexOf("q")!=-1)
                dPath = "M"+end+",2L"+start+",5L"+end+",8Z";
              else
                dPath = "M"+start+",2L"+end+",5L"+start+",8Z";
            }else{ 
              dPath = "M"+start+",2L"+end+",2L"+end+",8L"+start+",8Z";
            }
            return dPath;
          })
}

function svg2pdf(){

var pdfWidth = width + margin.left + margin.right,
    pdfHeight = height + margin.top + margin.bottom,
    doc = new jsPDF(pdfWidth>pdfHeight?"l":"p","pt",[pdfWidth, pdfHeight]),
    svgDoc = d3.select("svg").node();

doc.setDrawColor(0);
doc.setLineWidth(1);

d3.selectAll(".chromosome .bands>path:not(.acen)").each(function(){
	var coords = getCoords(this, svgDoc),
	  color = d3.rgb(d3.select(this).style("fill"));
	doc.setFillColor(color.r,color.g,color.b);
	doc.rect(coords.x, coords.y, coords.x2-coords.x, coords.y2-coords.y, 'FD');
});

d3.selectAll(".chromosome path.acen").each(function(d,i){
	var coords = getCoords(this, svgDoc),
	  color = d3.rgb(d3.select(this).style("fill"));
	doc.setFillColor(color.r,color.g,color.b);
	if(i%2==0)
	  doc.triangle(coords.x, coords.y, coords.x2, coords.y, (coords.x+coords.x2)/2, coords.y2, 'FD');
	else
	  doc.triangle(coords.x, coords.y2, coords.x2, coords.y2, (coords.x+coords.x2)/2, coords.y, 'FD');
});

doc.setFillColor(0,144,0);
d3.selectAll(".chromosome path.mapTrack").each(function(){
  var coords = getCoords(this, svgDoc);
  doc.rect(coords.x, coords.y, coords.x2-coords.x, coords.y2-coords.y, 'F');
});

d3.selectAll(".chromosome g.brush rect.extent").each(function(){
  if(d3.select(this).style("display")!="none"){
	var coords = getCoords(this, svgDoc);
	doc.setDrawColor(0,0,255);
	doc.rect(coords.x, coords.y, coords.x2-coords.x, coords.y2-coords.y, 'D');
  }
});

d3.selectAll(".track>rect:first-child").each(function(){
  var coords = getCoords(this, svgDoc),
      color = d3.rgb(d3.select(this).style("fill"));
  doc.setFillColor(color.r,color.g,color.b);
  doc.rect(coords.x, coords.y, coords.x2-coords.x, coords.y2-coords.y,'F');
})

doc.setDrawColor(0);

d3.selectAll(".track>.segment>path, .track.score>path, .track>.values>path").each(function(){
  var self = d3.select(this),
      coords = getCoords(this, svgDoc),
      opc = parseFloat(self.style("fill-opacity")),
      color = d3.rgb(self.style("fill")),
      d = self.attr("d").replace(/M|Z/g,"").split("L"),
      acc = [],
      d0 = d[0].split(/[,| ]/).filter(function(a){ return a.length>0; }),
      cx = +d0[0],
      cy = +d0[1];
  if(this.parentNode.parentNode.classList.value=="track exons")
    coords.y = coords.y + 6;
  for(var i=1; i<d.length; i++){
    var di = d[i].split(/[,| ]/).filter(function(a){ return a.length>0; }),
        dx = +di[0]-cx,
        dy = +di[1]-cy;
    acc.push([dx, dy]);
    cx += dx;
    cy += dy;
  }
  if(opc<1)
    color = applyOpacity(color,opc,{r:255,g:255,b:255});
  doc.setFillColor(color.r,color.g,color.b);
  doc.lines(acc, coords.x, coords.y, [1,1], 'F', true);
})

d3.selectAll(".track.vcfmaster>.segment>rect").each(function(){
  var coords = getCoords(this, svgDoc),
      color = d3.rgb(d3.select(this).style("fill"));
  doc.setFillColor(color.r,color.g,color.b);
  doc.roundedRect(coords.x, coords.y, coords.x2-coords.x, coords.y2-coords.y, 2, 2,'F')
})

doc.setFillColor(255);
doc.rect(0, 0, margin.left, height + margin.top + margin.bottom, 'F');
doc.rect(margin.left+width, 0, margin.left, height + margin.top + margin.bottom, 'F');

d3.selectAll("line").each(function(){
	var coords = getCoords(this, svgDoc);
	doc.line(coords.x, coords.y, coords.x2, coords.y2);
});

d3.selectAll(".x.axis path").each(function(){
	var coords = getCoords(this, svgDoc);
	doc.setDrawColor(0);
	doc.line(coords.x, coords.y2, coords.x2, coords.y2);
});

d3.selectAll("text").each(function(){
	var self = d3.select(this),
        coords = getCoords(this, svgDoc),
	t = d3.select(this.childNodes[0]).text(),
	color = d3.rgb(self.style("fill")),
	s = parseInt(self.style("font-size")),
	rotate = d3.transform(self.attr("transform")).rotate;
	doc.setTextColor(color.r,color.g,color.b);
	doc.setFontSize(s);
	doc.text(coords.x, coords.y+8, t, null, rotate?-rotate:null);
});

doc.save(d3.select("head>title").text()+".pdf");
}
