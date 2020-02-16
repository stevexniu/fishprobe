var body,
    php = typeof SQL == 'undefined',
    dataTracks = [],
    chromosomes,
    chr,
    start,
    end;

var fixedBounds = false;

var margin = 0, width = 0, height = 0;

if(php){
  if(window.location.protocol=="file:")
    d3.select("html")
      .append("div")
      .html("<p>This genome browser has been configurated to be hosted in a web server.</p><p>Please, set server parameter to false in the creation function to enable the local visualization mode.</p>")
      .style({"font-family":"sans-serif",margin:"0 auto","margin-top":"100px",padding:"10px","max-width":"900px"})
  else
    window.onload = function(){
      displayGenomeBrowser();
    }
}else{
  load_db(displayGenomeBrowser);
}

function displayGenomeBrowser(){
  body = d3.select("body");
  chromosomes = JSON.parse(d3.select("#chromosomes").text());
  execDB("SELECT trackid, trackname, type, color, data FROM tbl_tracks", function(res){
    dataTracks = res;
    d3.select("head").append("style")
      .text("body {font-family: sans-serif; font-size: 14px; color: #333;}"+
      "ul { list-style-type: none; padding: 0; }"+
      "pre { font-family: monospace; }"+
      "div.window { font-size: 12px; }"+
      "span.highlight { background-color: #ffa; }"+
      "div.tooltip { background-color: #fff; font-size: 12px; padding: 0px 14px; border-radius: 6px; border: 2px solid #ddd;}"+
      "input, button { border: 1px solid #aaa; border-radius: 3px; margin: 2px; padding: 2px 4px; font-size: 14px; color: #333; }"+
      "button, input[type=button] { background: #eee; background: linear-gradient(#fff, #eee); color: #555; cursor: pointer; }"+
      "button:disabled { color: #aaa; cursor: default; }")
    d3.select("html").on("click", function(){
        d3.select("div.tooltip").remove();
        d3.select("select").remove();
      })
    if(/^(\?r=)([\_\.0-9a-zA-Z]+)\:[0-9]+\-[0-9]+$/.test(window.location.search)){
      displayTracks(window.location.search.substring(3));
    } else {
      genomemap();
    }
  });
}

function execDB(sql, callback){
  if(php){
    d3.json('query.php')
      .header("Content-type", "application/x-www-form-urlencoded")
      .post("q="+encodeURI(sql), function(error,res){
        if (error) throw error;
          res = (res.length == 0)? false : res;
          callback(res);
      });
  }else{
    get_query(sql,callback);
  }
}

function gotoLocation(){
  var str = d3.select("input[name=r]").node().value,
      patt = new RegExp("^([\_\.0-9a-zA-Z]+)\:[0-9]+(\-[0-9]+)?$");
  if(patt.test(str)){
    displayTracks(str);
  }else{
    var option = d3.select("select[name=r]>option:first-child");
    if(!option.empty())
      displayTracks(option.property("value"));
      d3.select("select[name=r]").remove();
  }
  return false;
}

function searchGene(){
  var str = d3.select("input[name=r]").node().value,
  patt = new RegExp("^[\_\.0-9A-Za-z\-]{2,}$");
  d3.select("select[name=r]").remove();
  if(patt.test(str)){
    execDB("SELECT * FROM aux_genes WHERE name LIKE '"+str+"%' LIMIT 100", function(gene){
      if(gene != false){
        var BCR = d3.select("input[name=r]").node().getBoundingClientRect();
        body.append("select")
              .attr("name","r")
              .attr("size", gene.length+1)
              .on("change", function(){
                displayTracks(this.value);
                d3.select(this).remove();
              })
              .style({"position":"absolute", "top":(BCR.bottom+window.scrollY)+"px", "left": (BCR.left)+"px"})
              .selectAll("option").data(gene)
              .enter().append("option")
                .property("value", function(d){ return d[0]+":"+d[1]+"-"+d[2]; })
                .text(function(d){ return (d[4])?(d[3]+" - "+d[4].replace(/\<[^<]+\>/g,"")):d[3]; });
      }
    });
  }
}

var anchors = {
  uniprot: 'http://www.uniprot.org/uniprot/',
  interpro: 'https://www.ebi.ac.uk/interpro/entry/',
  pfam: 'http://pfam.xfam.org/family/',
  symbol: 'https://www.ncbi.nlm.nih.gov/gene/?term=',
  geneid: 'https://www.ncbi.nlm.nih.gov/gene/',
  ensembl: 'http://www.ensembl.org/id/',
  go: 'https://www.ebi.ac.uk/QuickGO/GTerm?id=',
  keggs: 'http://www.genome.jp/dbget-bin/www_bget?',
  ez: 'http://www.brenda-enzymes.org/enzyme.php?ecno=',
  biosystems: "https://www.ncbi.nlm.nih.gov/biosystems/?term=",
  ips: 'https://www.ebi.ac.uk/interpro/signature/'
}

var replaceAnchors = function(d){
  if(d.match(/^[a-z]+\:/)){
    for(term in anchors){
      if(d.match(new RegExp("^"+term+"\:")))
        return d.replace(new RegExp("^"+term+":[A-Za-z0-9\.\_\:]+"),function(o){ 
          o = o.replace(term+":","");
          return '<a target="_blank" href="'+anchors[term]+o+'">'+o+'</a>';
        });
    }
  }
  return d;
}

function searchGlobal(){
  var str = d3.select("input[name=r]").node().value,
      patt = new RegExp("^[\_\.0-9A-Za-z\-]{2,}$"),
      highlighter = function(s){
        s = replaceAnchors(s);
        return s.replace(new RegExp("(?![^<]*>)"+str, 'gi'),function(o){ return '<span class="highlight">'+o+'</span>'; });
      }

  if(patt.test(str)){
    execDB("SELECT * FROM aux_genes WHERE name LIKE '%"+str+"%' OR score LIKE '%"+str+"%' LIMIT 100", function(res){
      if(res != false){
        var win = showWindow(),
            li = win.append("ul")
                  .selectAll("li").data(res)
                .enter().append("li")
        li.append("h5")
          .html(function(d){ return highlighter(d[3]); })
          .style({"cursor":"pointer","color":"#1f77b4","display":"inline"})
          .on("click", function(d){
            chr = d[0];
            start = d[1];
            end = d[2];
            displayTracks();
          })
        li.each(function(d){
          var dat = d[4]?d[4].split("|"):[];
          if(dat.length > 1){
            d3.select(this).append("ul")
              .style('padding','2px 0 8px 10px')
              .selectAll("li")
              .data(dat)
            .enter().append("li")
              .html(function(d){ return highlighter(d); })
          }
          if(dat.length == 1){
            d3.select(this).append("span")
              .html(" - " + highlighter(dat[0]))
          }
        })
      }else{
            showWindow("not found");
      }
    });
  }
}

function imgInfo(sel, callback){
    sel.append("img")
	.attr("alt", "info")
	.attr("width", 14)
	.attr("height", 14)
	.attr("src", infoIcon)
	.style({"cursor": "pointer", "padding-left":"4px"})
	.on("click", function(){ callback(showWindow()); });
}

function displayTooltip(infoButton, displayInfo){
  if(d3.select("div.tooltip").empty()){
    d3.event.stopPropagation();
    var BCR = infoButton.getBoundingClientRect();
    body.append("div")
	.attr("class", "tooltip")
	.style({"position":"absolute",
	"top":(BCR.bottom+window.scrollY)+"px",
	"left": (BCR.right)+"px"})
	.call(displayInfo);
  }
}

function displaySelector(gChr, chrScale, selector){

  var brush = d3.svg.brush()
      .x(chrScale)
      .extent(selector?selector:[0,0])
      .on("brush", brushmove);

  var gBrush = gChr.select("g.brush")
  if(gBrush.empty()){

    gBrush = gChr.append("g")
      .attr("class","brush")
      .call(brush);

    gBrush.select(".background")
      .attr("height", 10)
      .style("cursor","pointer")
      .on("mousedown", !selector ? function(){
        var id = gChr.attr("id");
        d3.selectAll(".brush .extent, .brush .resize").style("display",function(){ return id==this.parentNode.parentNode.id?null:"none"; });
      } : null)

    gBrush.selectAll(".resize")
      .style("cursor","move")
      .append("rect")
        .attr("class", "handle")
        .attr("height", 10)
        .attr("width", 2)
        .style("fill","blue")

    gBrush.select(".extent")
      .attr("height", 10)
      .style({"fill-opacity": ".8", "stroke": "blue", "fill": "none"});

  }else{

    gBrush.call(brush);

  }

  function brushmove() {
    var extent = brush.extent(),
        bp = (typeof start != 'undefined')?end-start:1e+5;

    if(bp>chrScale.domain()[1]){
      start = 0;
      end = chrScale.domain()[1];
      gBrush.call(brush.extent([start,end]));
    }else{
      start = Math.round(d3.min(extent));
      end = Math.round(d3.max(extent));
      if(end-start==0){
        end = start+(bp/2);
        start = start-(bp/2);
        if(start<0){
          start = 0;
          end = bp;
        }
        if(end>chrScale.domain()[1]){
          start = chrScale.domain()[1]-bp;
          end = chrScale.domain()[1];
        }
        gBrush.call(brush.extent([start,end]));
      }
    }
      
    if(typeof selector != 'object'){
      body.select("a>button").attr("disabled", null)
        .on("click", function(){
          chr = gChr.attr("id");
          displayTracks();
        })
    }else{
      displayTracks();
    }
  }
}

function getTrackY(svg){
  try{
    var bbox = svg.node().getBBox(),
    trackY = bbox.y+bbox.height+20;
    return trackY;
  }catch(e){
    return -1;
  }
}

function resizeHeight(sel){
  var bbox = sel.node().getBBox();
  height =  Math.ceil(bbox.height);
  d3.select("svg").attr("height", height + margin.top + margin.bottom);
  d3.select("#mask rect").attr("height", height);
}

function addTitle(sel){
  sel.append("title")
    .text(function(d){
      var txt = "";
      if(d[5])
        txt += d[5]+"\n\n";
      if(d[6])
        txt += d[6].replace(/\|/g,"\n").replace(/\<[^<]+\>/g,"")+"\n\n";
      txt += d[2]+":"+(d[3]+1)+(((d[3]+1)==d[4])?"":("-"+d[4]));
      return txt;
    });
}

function showVCFtext(d,data){
  data = data.info?data.info:data;

  var sel = showWindow();
  sel.append("style").text('div#vcfdata span.attr { color: #888; }');
  sel = sel.append("div").attr("id",'vcfdata');
  
  if(d[5])
    sel.append('h4').text(d[5])
  if(d[6]){
    sel.append('ul').selectAll('li')
        .data(d[6].split("|"))
      .enter().append('li')
        .each(function(d){
          d = d.split('=');
          var self = d3.select(this);
          self.append('span')
            .attr('class','attr')
            .text(d[0] + ': ')
            .attr('title',data[d[0]])
          self.append('span')
            .text(d.length==1?d[0]:d[1])
        })
  }
  sel.append('p').text(d[2]+":"+(d[3]+1)+(((d[3]+1)==d[4])?"":("-"+d[4])))
}

function keypressed() {
  if(d3.select("input:focus").empty()){
    if(d3.event.code=="ArrowLeft"){
          var move = parseInt((end-start)*0.1);
          start = start - move;
          end = end - move;
          displayTracks();
    }
    if(d3.event.code=="ArrowRight"){
          var move = parseInt((end-start)*0.1);
          start = start + move;
          end = end + move;
          displayTracks();
    }
  }
}

function downloadFasta(q){ 
  var pre = document.querySelectorAll('div.window pre'),
      texto = pre[q].textContent,
      name = texto.substring(1,texto.indexOf("\n"));
  var blob = new Blob([texto], {type: 'text/plain'});
  fileDownload(blob,name+".fa");
}

function changeFasta(){
  var pre = document.querySelectorAll('div.window pre'),
      but = document.querySelectorAll('div.window button');
  if(pre[0].style.display == 'none'){
    but[0].onclick = function() { downloadFasta(0) };
    but[2].onclick = function() { blast(0) };
    pre[0].style.display = 'block';
    pre[1].style.display = 'none';
  }else{
    but[0].onclick = function() { downloadFasta(1) };
    but[2].onclick = function() { blast(1) };
    pre[1].style.display = 'block';
    pre[0].style.display = 'none';
  }
}

function blast(q) {
  var pre = document.querySelectorAll('div.window pre'),
      ncbiBlastLink = 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&PROGRAM=blast'+(q?'p':'x')+'&DATABASE=nr&QUERY='+encodeURI(pre[q].textContent);
      var link = document.createElement('a');
      link.href = ncbiBlastLink;
      link.target = '_blank';
      var clicEvent = new MouseEvent('click', {
        'view': window,
        'bubbles': true,
        'cancelable': true
      });
      link.dispatchEvent(clicEvent);
}

var xhrGetTracks;
function getTracks(chr, start, end, callback){
  if(php){
    if(typeof xhrGetTracks != 'undefined')
      xhrGetTracks.abort();
    xhrGetTracks = d3.json('query.php')
	.header("Content-type", "application/x-www-form-urlencoded")
	.post("r="+encodeURI(chr+":"+start+"-"+end), function(error,res){
            if (error) throw error;
	    callback(res);
        });
  }else{
    callback(get_tracks(chr, start, end));
  }
}

function genomebrowser(){

  margin = {top: 80, right: 50, bottom: 0, left: 50};
  width = 1150 - margin.left - margin.right;
  height = 800 - margin.top - margin.bottom;

  body.html("")
    .on("keydown", keypressed);

  var nav = body.append("div")
    .style({"width":(width + margin.left + margin.right)+"px", "margin": "0 auto"})
    .append("div")
      .style({"padding-left": margin.left+"px", "padding-right": margin.right+"px", "margin-top": "20px"})
  nav.append("span")
       .attr("class","location")
  nav.append("input")
     .attr("type", "text")
     .attr("name", "r")
     .style("width","200px")
     .on("keyup", function(){
        if (d3.event.keyCode == 13) { gotoLocation(); }
	else { searchGene(); }
      })
     .on("paste", function(){
	searchGene();
      })
  nav.append("button")
     .text("Go to")
     .on("click", gotoLocation);
  nav.append("button")
     .text("Search")
     .on("click", searchGlobal);

  imgInfo(nav, function(sel){
    execDB("SELECT * FROM aux_genes WHERE _ROWID_ = (abs(random()) % (SELECT max(_ROWID_) FROM aux_genes)) LIMIT 1", function(gene) {
      var example = gene[0];
      sel.append("style").text(".window p { margin: 20px 0; } .window p>svg { vertical-align: middle; }")
      sel.append("p").html("Input a region specified like '"+(example[0]+":"+example[1]+"-"+example[2])+"' or '"+(example[3])+"', then click <button>Go to</button>")
      sel.append("p").html("Input a word to search in browser, then click <button>Search</button>")
      sel.append("p").html(infoImages.clickHide+" Click on a track's name to hide/show it.")
      sel.append("p").html("Drag n' drop (horizontal) or right and left arrows to advance in the browser.")
      sel.append("p").html("Drag n' drop (vertical) to reorder tracks.<span style=\"padding: 0 4%;\"></span>"+infoImages.move)
      sel.append("p").html(infoImages.rescale+"Click on a value track's scale to change it.")
      sel.append("p").html("Click on genes, exons, domains or vcf will display a window with information about the feature.")
      sel.append("p").html("Hover the features for information.<span style=\"padding: 0 4%;\"></span>"+infoImages.clickInfo)
      sel.append("p").html("In vcf's information window, hover the acronyms for explanation.")
    });
  });

  if(dataTracks.map(function(d){ return d[2]; }).indexOf("value")!=-1){
    nav.append("button")
      .attr("class","bounds")
      .text("BOUNDS")
      .attr("title","set fixed boundaries for 'value' tracks")
      .style("float","right")
      .on("click", function(){
        fixedBounds = !fixedBounds;
        d3.select(this).style("background",fixedBounds?"#aaa":null)
        displayTracks();
      });
  }

  nav.append("button")
    .text("MAP")
    .attr("title","go back to genomemap")
    .style("float","right")
    .on("click", function(){
      window.location.href = window.location.href.split("?")[0];
    });

  nav.append("button")
    .text("FASTA")
    .attr("title","get FASTA of this region")
    .style("float","right")
    .on("click", function(){
          var queryCallback = function(error, data) {
            if (error) throw error;
            if(data && data.dna){
              var win = showWindow();
              win.append("button").text("download").attr("onclick","downloadFasta(0);");
              win.append("pre").text(">" + chr + ":" + start + "-" + end + "\n" + data.dna);
            }else
              showWindow("No sequence available");
          }
          if(php)
	    d3.json("./query.php?r="+chr+":"+start+"-"+end, queryCallback);
          else
            get_dna([0,0,chr,start,end], queryCallback);
    });

var svg = body.select("body>div").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom);

  svg.append("style")
    .text("text { font: 10px sans-serif; } "+
".axis text, text.trackName, text.noData, text.legend { fill: #333; }"+
"text.trackName { cursor: pointer; }"+
".segment text { font: 8px monospace; } "+
".track.vcf .segment text { fill: #fff; }"+
".axis path, .axis line { fill: none; stroke: #333; shape-rendering: crispEdges; } "+
".buttons {fill: grey;}" +
".buttons>circle { cursor: pointer; opacity: 0.2; }"+
".bands>* {stroke: #000;  fill: #FFF;} "+
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

  svg.append("defs")
    .append("clipPath")
    .attr("id", "mask")
    .append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", width)
        .attr("height", height);

  var buttons = svg.append("g")
      .attr("class", "buttons")
      .attr("transform", "translate("+margin.left+","+(margin.top-64)+")")

  buttons.append("polygon")
      .attr("points", "0,5 5,0 5,5 10,0 10,10 5,5 5,10")
  buttons.append("circle")
      .attr("cx","6")
      .attr("cy","5")
      .attr("r","7")
      .on("click", function(){
          var move = end-start;
          start = start-move;
          end = end-move;
          displayTracks();
        })

  buttons.append("polygon")
      .attr("points", "0,5 10,0 10,10")
      .attr("transform", "translate("+15+",0)")
  buttons.append("circle")
      .attr("transform", "translate("+15+",0)")
      .attr("cx","6")
      .attr("cy","5")
      .attr("r","7")
      .on("click", function(){
          var move = parseInt((end-start)*0.1);
          start = start-move;
          end = end-move;
          displayTracks();
        })

  buttons.append("rect")
      .attr("transform", "translate("+(width/2-20)+",0)")
      .attr("x","0")
      .attr("y","3.5")
      .attr("width","10")
      .attr("height","3")
  buttons.append("circle")
      .attr("transform", "translate("+(width/2-20)+",0)")
      .attr("cx","5")
      .attr("cy","5")
      .attr("r","7")
      .on("click", function(){
          var bp = end-start,
              center = start+parseInt(bp/2),
              zoomout = parseInt(bp/0.8);
          start = parseInt(center-zoomout/2);
          end = parseInt(center+zoomout/2);
          displayTracks();
        })

  buttons.append("path")
      .attr("d", "M 3.6308594 0 A 3.6302359 3.6302359 0 0 0 0 3.6308594 A 3.6302359 3.6302359 0 0 0 3.6308594 7.2597656 A 3.6302359 3.6302359 0 0 0 5.4355469 6.7753906 L 8.6601562 10 L 10 8.6601562 L 6.7753906 5.4355469 A 3.6302359 3.6302359 0 0 0 7.2597656 3.6308594 A 3.6302359 3.6302359 0 0 0 3.6308594 0 z M 3.6308594 0.90820312 A 2.6832178 2.6832178 0 0 1 6.3125 3.5917969 A 2.6832178 2.6832178 0 0 1 3.6308594 6.2753906 A 2.6832178 2.6832178 0 0 1 0.94726562 3.5917969 A 2.6832178 2.6832178 0 0 1 3.6308594 0.90820312 z ")
      .attr("transform", "translate("+(width/2-5)+",0)")

  buttons.append("path")
      .attr("d", "M 3.5 0 L 3.5 3.5 L 0 3.5 L 0 6.5 L 3.5 6.5 L 3.5 10 L 6.5 10 L 6.5 6.5 L 10 6.5 L 10 3.5 L 6.5 3.5 L 6.5 0 L 3.5 0 z ")
      .attr("transform", "translate("+(width/2+10)+",0)")
  buttons.append("circle")
      .attr("transform", "translate("+(width/2+10)+",0)")
      .attr("cx","5")
      .attr("cy","5")
      .attr("r","7")
      .on("click", function(){
          var bp = end-start,
              center = start+parseInt(bp/2),
              zoomin = parseInt(bp*0.8);
	  if(zoomin<100)
	    zoomin = 100;
          start = parseInt(center-zoomin/2);
          end = parseInt(center+zoomin/2);
          displayTracks();
        })

  buttons.append("polygon")
      .attr("points", "0,5 -10,0 -10,10")
      .attr("transform", "translate("+(width-15)+",0)")
  buttons.append("circle")
      .attr("transform", "translate("+(width-15)+",0)")
      .attr("cx","-6")
      .attr("cy","5")
      .attr("r","7")
      .on("click", function(){
          var move = parseInt((end-start)*0.1);
          start = start+move;
          end = end+move;
          displayTracks();
        })

 buttons.append("polygon")
      .attr("points", "0,5 -5,0 -5,5 -10,0 -10,10 -5,5 -5,10")
      .attr("transform", "translate("+width+",0)")
  buttons.append("circle")
      .attr("transform", "translate("+width+",0)")
      .attr("cx","-6")
      .attr("cy","5")
      .attr("r","7")
      .on("click", function(){
          var move = end-start;
          start = start+move;
          end = end+move;
          displayTracks();
        })

  svg = svg.append("g")
      .attr("class", "tracks")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  svg.append("g")
      .attr("class","chromosome")
      .attr("transform", "translate(0," + (-40) + ")");

  var tracks = svg.append("g")
      .attr("clip-path", "url(#mask)");

  svg.append("g")
      .attr("class", "x axis");

  var track = tracks.selectAll(".track")
      .data(dataTracks)
    .enter().append("g")
      .attr("class", function(d){ return "track "+d[2]; })
      .style("fill", function(d){ return d[3]; })
      .property("scale", function(d){ return d[2]=="score"?JSON.parse(d[4]):null; })

  track.append("rect")
          .attr("class","rectBack")
          .attr("x",0)
          .attr("y",-20)
          .attr("width",width)
          .attr("height",0)

  track.append("text")
        .attr("class","trackName")
        .text(function(d){ return d[1]; })
        .attr("y", -6)
        .attr("x", 6)
        .on("click",function(d){
          d[5] = d[5]?0:1;
          displayTracks();
        })

  var initX;
  var widthScale = d3.scale.linear()
        .domain([0,width]);
  var sortTracks = function(a,b){
    a = a[6];
    b = b[6];
    return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
  }

  var drag = d3.behavior.drag()
             .on('dragstart', function(){
               tracks.style("cursor","move");
               initX = d3.mouse(tracks.node())[0];
             })
             .on('drag', function(thisData) {
               var coor = d3.mouse(tracks.node())[1];
               if(coor%10==0){
                 track.each(function(d){
                     if(coor > d[6] && coor < d[6]+d[7]){
                       var idY = thisData[6];
                       thisData[6] = d[6];
                       d[6] = idY;
                       track.sort(sortTracks)
                         .style("display","none")
                         .each(function(d,i){
                           d[6] = getTrackY(tracks);
                           d3.select(this)
                               .attr("transform","translate(0,"+d[6]+")")
                               .style("display",null)
                               .select(".rectBack")
                                 .style("fill",i%2==1?"#f0f0f0":"#fcfcfc")
                         })
                     }
                 });
               }
             })
             .on('dragend', function() {
               tracks.style("cursor",null);
               var coor = d3.mouse(tracks.node())[0],
                   x = coor - initX;
               if(Math.abs(x)>20){
                 widthScale.range([0,end-start]);
                 var move = parseInt(widthScale(x));
                 start = start - move;
                 end = end - move;
                 displayTracks();
               }
             });

  track.call(drag)

  displayButtons();
  bioinfoLogo();

  return svg;
}

function decodeLocation(strLocation){
  chr = strLocation.split(":")[0];
  start = parseInt(strLocation.split(":")[1].split("-")[0]);
  end = parseInt(strLocation.split(":")[1].split("-")[1]);

  if(isNaN(end)){
    if(start>50){
      start = start-50;
      end = start+100;
    }else{
      start = 0;
      end = 100;
    }
  }
}

function displayTracks(strLocation){

  if(strLocation)
    decodeLocation(strLocation);

  var svg = d3.select("svg g.tracks");
  if(svg.empty())
    svg = genomebrowser();

  var bp = end-start;
  if(bp<100){
    end = start+100;
    bp = 100;
  }

  if(start<0){
    start = 0;
    end = bp;
  }

  var json = chromosomes[chr.toLowerCase()];

  var chrEnd = d3.max(json, function(d){
      return d[1];
  });

  if(end>chrEnd){
    end = chrEnd;
    start = chrEnd-bp;
    if(start<0){
      start = 0;
      bp = end;
    }
  }

var chrScale = d3.scale.linear()
    .range([0, width])
    .domain([0, chrEnd]);

var x = d3.scale.linear()
    .range([0, width])
    .domain([start, end]);

var scoreScale = d3.scale.linear()
    .clamp(true);

var valueScale = d3.scale.linear()
    .range([50, 0])
    .clamp(true);

var displayValueTrack = function(track){
  x.domain([start,end]);
  valueScale.domain(track.property("scale")).nice();
  var vAxis = d3.svg.axis()
	.scale(valueScale)
	.tickValues(valueScale.domain())
	.orient("left");
  track.selectAll("g.values>path")
        .attr("d", function(d){
          var s = x(d[3]),
              e = x(d[4]),
              h = valueScale(d[6]),
              b = valueScale(0);
          if(e-s<1){
            s = ((s+e)/2)-.5;
            e = ((s+e)/2)+.5;
          }
          return "M"+s+","+h+"L"+e+","+h+"L"+e+","+b+"L"+s+","+b+"Z"; }
        );
  track.select(".axis").call(vAxis);
}

var cutSegmentStart = function(x){
  var bound = -margin.left;
  return ((x<bound)?bound:x)
}

var cutSegmentEnd = function(x){
  var bound = width+margin.right;
  return ((x>bound)?bound:x)
}

var xAxis = d3.svg.axis()
    .scale(x)
    .orient("top")
    .tickFormat(function(d){
      if(bp>6666){
        if(d>=1e+6)
          return d3.round((d/1e+6),3)+"M";
        if(d>=1e+3)
          return d3.round((d/1e+3),3)+"K";
      }
      return d3.format(",")(d);
    });

  var bpFormat = function(bp){
    if(bp%1e+6==0)
      return bp/1e+6 + "Mb";
    if(bp%1e+3==0)
      return bp/1e+3 + "kb";
    return bp + "bp";
  }

  d3.select("span.location")
      .text("location "+chr+":"+start+"-"+end+" "+bpFormat(bp)+" ")

  d3.select("input[name=r]")
      .property("value","");

  if(typeof displayChr == 'function'){
    var gChr = svg.select("g.chromosome").attr("id", chr);
    displayChr(gChr, chrScale, json);
    displaySelector(gChr, chrScale, [start, end]);
  }

  svg.select(".x.axis")
      .call(xAxis);

  var gTracks = svg.select("g[clip-path]")

  var gTrack = gTracks.selectAll(".track")
        .style("display","none");

  getTracks(chr,start,end,function(allData){

    gTrack.each(function(d,trackIndex){

      var trackID = d[0], 
      type = d[2],
      hiddenTrack = d[5],
      trackY = d[6] = getTrackY(gTracks);

      var data = allData['track'+d[0]];

      var track = d3.select(this)
        .attr("transform", "translate(0,"+trackY+")")
        .style("display",null);

      if(!data||data[0][2]==false||hiddenTrack){
        track.selectAll(":not(.trackName):not(.rectBack):not(.noData)").remove();
        var noData = track.selectAll(".noData")
          .data([hiddenTrack?"(hidden track)":(data?"(zoom in to view this track)":"(empty region)")])
        noData.enter().append("text")
          .attr("class","noData")
          .attr("transform", "translate("+(width/2)+",0)")
          .attr("text-anchor", "middle")
        noData.text(String);

      }else{

      track.select(".noData").remove();

      switch(type){
        case "value":

          if(fixedBounds){
            if(this.newscale)
              this.scale = this.newscale;
            else
              this.scale = JSON.parse(d[4]);
          }else{
            delete this.newscale;
            var minval = d3.min(data,function(d){ return +d[6]; }),
                maxval = d3.max(data,function(d){ return +d[6]; });
            if(minval>0)
              minval = 0;
            if(maxval<0)
              maxval = 0;
            this.scale = [minval,maxval];
          }

          track.selectAll(".track>g")
            .data(["values","axis"])
            .enter().append("g")
            .attr("class",String)
            .each(function(d){
              if(d=="axis"){
                d3.select(this).attr("transform", "translate("+(width-1)+")")
                .on("click", function(){
		  if(d3.select("div.tooltip").empty()){
		    d3.event.stopPropagation();
		    var BCR = this.getBoundingClientRect(),
		    tooltip = body.append("div")
			.attr("class", "tooltip")
			.style({"position":"absolute",
			"top":(BCR.top+window.scrollY)+"px",
			"left":(BCR.right)+"px"})
			.on("click", function(){ d3.event.stopPropagation(); })
		    tooltip.append("p")
			.append("input")
			  .attr("type","number")
			  .attr("name","max")
			  .style("width", "50px")
			  .property("value",this.childNodes[1].textContent.replace(",",""))
		    tooltip.append("p")
			.append("input")
			  .attr("type","number")
			  .attr("name","min")
			  .style("width", "50px")
			  .property("value",this.childNodes[0].textContent.replace(",",""))
		    tooltip.append("p")
			.append("input")
			  .attr("type", "button")
			  .attr("name", "setScale")
			  .style("width", "50px")
			  .property("value", "set")
			  .on("click", function(){
			    fixedBounds = true;
			    d3.select("button.bounds").style("background","#aaa")
			    var newScale = [d3.select("input[name=min]").property("value"),d3.select("input[name=max]").property("value")];
			    track.property("newscale",newScale);
			    track.property("scale",newScale);
			    displayValueTrack(track);
			  })
		  }
                });
              }
            });

          var values = track.select(".values").selectAll("path")
                .data(data, function(d){ return d[0]; })
          values.enter().append("path").call(addTitle);
          values.exit().remove();

          displayValueTrack(track);

        break;
        case "score":

          scoreScale.domain(track.property("scale")).nice();

          var scores = track.selectAll("path")
              .data(data, function(d){ return d[0]; })
          scores.enter().append("path")
            .style("fill-opacity",function(d){ return scoreScale(d[6]); })
            .call(addTitle);
          scores.exit().remove();
          scores.attr("d",function(d){ 
                return "M"+x(d[3])+",0L"+x(d[4])+",0L"+x(d[4])+",10L"+x(d[3])+",10Z";
              })

        break;
        case "vcf":

          var info = JSON.parse(d[4]),
              show = info.show ? info.show : [],
              dataLegend = ["Variants"],
              infoColor = d3.scale.ordinal()
                .range(["#ff7f0e","#2ca02c","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"]);

          show.forEach(function(d){
            dataLegend.push(d);
          })

          track.selectAll(".legend")
            .data(dataLegend)
              .enter().append("text")
            .attr("class","legend")
            .attr("x",6)
            .attr("y",function(d,i){ return 10+(12*i); })
            .style("font-size", "8px")
            .text(String)

          var segment = track.selectAll(".segment")
            .data(data, function(d){ return d[0]; })

          var smgEnter = segment.enter().append("g")
            .attr("class", "segment")

          segment.exit().remove();

          smgEnter.each(function(d){
            var thiz = d3.select(this),
                attrs = d[6].split("|"),
                aux = [],
                score = {};

            attrs.forEach(function(p){
              aux = p.split("=");
              if(aux.length==1)
                score[p] = p;
              else
                score[aux[0]] = aux[1];
            })

            thiz.append("path")
              .attr("d","M-1,0L2,0L2,10L-1,10Z")
              .style("fill","steelblue")

            if(score['AF']!=null){
              thiz.append("path")
                .attr("d",function(d){
                  var h = 10*(1-score['AF']);
                  return "M-1,"+h+"L2,"+h+"L2,10L-1,10Z";
                })
                .style("fill","firebrick")
            }

            var h = 22,
                val = "";
            show.forEach(function(name){
                  val = score[name];
                  if(val){
                    thiz.append('rect')
                      .attr("x",-3)
                      .attr("height",10)
                      .attr("y",h-8)
                      .attr("width",7)
                      .attr('rx',2)
                      .style('fill',function(){
                        var patt = new RegExp("_MAF$");
                        if(patt.test(name) || name == "GMAF"){
                          var n = parseFloat(val.split(":")[1])
                          return d3.lab((1-n)*100,0,-40).toString();
                        }
                        if(name == "IMPACT"){
                          return val=="HIGH"?"#a50f15":(val=="MODERATE"?"#de2d26":(val=="LOW"?"#fb6a4a":"#fc9272"));
                        }
                        return infoColor(name);
                      })
                    thiz.append('text')
                      .attr("x",-2)
                      .attr("y",h)
                      .text(val.slice(0,1).toUpperCase())
                      .append("title")
                        .text(name+'='+val)
                  }
                  h = h + 12;
            })
          })

          smgEnter.style('cursor','pointer')
            .on("click",function(d){ showVCFtext(d,info); });

          segment.attr("transform",function(d){ return "translate("+x(d[3])+",0)"; })

        break;
        case "vcfsample":

          var segment = track.selectAll(".segment")
            .data(data, function(d){ return d[0]; })

          var smgEnter = segment.enter().append("g")
            .attr("class", "segment")

          segment.exit().remove();

          smgEnter.append("path")
              .attr("d","M-0.5,0L2,0L2,10L-0.5,10Z")
              .style("fill","steelblue")

          smgEnter.append("path")
              .attr("d",function(d){
                var h = 10,
                  attr = d[6].split("|"),
                  AD = attr.filter(function(d){ return d.indexOf("AD=")==0; }),
                  RO = 1,
                  AO = 0;
                if(AD.length==1){
                  AD = AD[0].split("=")[1].split(",");
                  RO = parseInt(AD[0]);
                  AO = parseInt(AD[1]);
                }else{
                  AD = [];
                  AD[0] = attr.filter(function(d){ return d.indexOf("RO=")==0; }),
                  AD[1] = attr.filter(function(d){ return d.indexOf("AO=")==0; });
                  if(AD[0].length==1 && AD[1].length==1){
                    RO = parseInt(AD[0][0].split('=')[1]);
                    AO = parseInt(AD[1][0].split('=')[1]);
                  }
                }
                h = h-(h/(RO+AO)*AO);
                return "M-0.5,"+h+"L2,"+h+"L2,10L-0.5,10Z";
              })
              .style("fill","firebrick")

          var format = JSON.parse(d[4]);
          smgEnter.style('cursor','pointer')
            .on("click", function(d){ showVCFtext(d,format); })
            .call(addTitle);

          segment.attr("transform",function(d){ return "translate("+x(d[3])+",0)"; })

        break;
        default:

          var segment = track.selectAll(".segment")
            .data(data, function(d){ return d[0]; })

          var smgEnter = segment.enter().append("g")
            .attr("class", "segment")

          segment.exit().remove();

          smgEnter.append("path");
	  smgEnter.append("text")
	      .attr("y", 20)
	      .text(function(d){
		     var name = d[5]?d[5]:"";
		     if (d[7]=="+")
		       return name+">";
		     if(d[7]=="-")
		       return "<"+name;
		     return name;
		   });
	  smgEnter.call(addTitle);

          segment.selectAll("text")
	    .attr("x", function(d){ return x(d[3]); })

          segment.selectAll("path")
            .attr("d", function(d){
	      var sgmStart = cutSegmentStart(x(d[3])), sgmEnd = cutSegmentEnd(x(d[4]));

	      if(type == "domain")
	        return "M"+sgmStart+",0L"+(sgmStart+1)+",0L"+(sgmStart+1)+",4L"+(sgmEnd-1)+",4L"+(sgmEnd-1)+",0L"+sgmEnd+",0L"+sgmEnd+",10L"+(sgmEnd-1)+",10L"+(sgmEnd-1)+",6L"+(sgmStart+1)+",6L"+(sgmStart+1)+",10L"+sgmStart+",10Z";

	      if(type == "exons" && d[11]){
                var i;
                for(i = 12; i<=13; i++)
                  if(d[i].slice(-1)==",")
                    d[i] = d[i].slice(0,-1);
	        var dAttr = "M"+sgmStart+",6L"+sgmStart+",4L",
                    blockSizes = d[12].split(","),
                    blockStarts = d[13].split(",");
	        for(i = 0; i<d[11]; i++){
	          var exonStart = x(d[3]+parseInt(blockStarts[i])), exonEnd = x(d[3]+parseInt(blockStarts[i])+parseInt(blockSizes[i]));
	          if(exonEnd>0 && exonStart<width){
	            dAttr += exonStart + ",4L"+exonStart+",0L"+exonEnd+",0L"+exonEnd+",10L"+exonStart+",10L"+exonStart+",4L";
	          }
	        }
	        dAttr += sgmEnd+",4L"+sgmEnd+",6Z";
	        return dAttr;
	      }

	      return "M"+sgmStart+",0L"+sgmEnd+",0L"+sgmEnd+",10L"+sgmStart+",10Z";
            })

	    smgEnter.style('cursor','pointer').on("click", function(d){
                var queryCallback = function(error, data) {
                  if (error) throw error;
                  var win = showWindow();
                  win.append("h3").text(d[5])
                  win.append("span").text(d[2]+":"+d[3]+"-"+d[4])
                  win.append("br")
                  if(d[6])
                    win.append("ul").selectAll("li")
                        .data(d[6].split("|"))
                      .enter().append("li")
                        .html(replaceAnchors)
                  if(data && data.dna){
                    win.append("button").text("download").attr("onclick","downloadFasta(0);")
                    if(data.prot)
                      win.append("button").text("dna/protein").attr("onclick","changeFasta();")
                    win.append("button").text("blast").attr("onclick","blast(0);")
                    win.append("pre").text(">" + d[5] + "\n" + data.dna)
                    if(data.prot)
                      win.append("pre").style("display","none").text(">" + d[5] + "\n"+ data.prot);
                  }else
                    win.append("p").text("No sequence available");
                }
                if(php)
	          d3.json("./query.php?r="+d[2]+":"+d[3]+"-"+d[4]+"&name="+d[5]+((d[7]=="-")?"&strand=reverse":"")+((d[8]||d[9])?"&thickStart="+d[8]+"&thickEnd="+d[9]:"")+(d[11]?"&blockCount="+d[11]+"&blockSizes="+d[12]+"&blockStarts="+d[13]:""), queryCallback);
                else
                  get_dna(d, queryCallback);
	    })

          var lines = [];
          lines.push(-Infinity);
          segment.attr("transform",function(d){
            var selfCoords = this.getBBox(),
            i = 0;
            while(lines[i]>=selfCoords.x){
              i++;
              if(typeof lines[i] == 'undefined')
                lines.push(-Infinity);
            }
            lines[i] = selfCoords.x+selfCoords.width;
            return "translate(0,"+(i*24)+")";
          });
      }
      }
      var rectBack = track.select(".rectBack")
            .style("fill",trackIndex%2==1?"#f0f0f0":"#fcfcfc")
            .attr("height",0),
          trackBBox = track.node().getBBox(),
          trackHeight = d[7] = trackBBox.y+trackBBox.height+30;
      rectBack.attr("height",trackHeight);
    });
    resizeHeight(svg);
  });
}
