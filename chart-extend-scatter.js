
var lastEvent;
Chart.defaults.scatterWithBox = Chart.defaults.scatter;
var isDown=0;
var lastEvent2;
var e0={},e1={};
var originalDraw = Chart.controllers.scatter.prototype.draw
// var custom = Chart.controllers.scatter.extend({
var custom;
Chart.helpers.extend(Chart.controllers.scatter.prototype, {
    draw: function() {
        // Call super method first
        // Chart.controllers.scatter.prototype.draw.call(this, ease);
        // Chart.controllers.scatter.prototype.draw.apply(this,arguments);
        originalDraw.apply(this,arguments);
        // Chart.controllers.scatter.prototype.draw.apply(this,arguments);
        // call(this, ease);
        // Now we can do some custom drawing for this dataset. Here we'll draw a red box around the first point in each dataset
        var meta = this.getMeta();
        var pt0 = meta.data[0];
        var radius = pt0._view.radius;
        var ctx = this.chart.chart.ctx;

        ctx.save();
        ctx.strokeStyle = 'red';
        ctx.lineWidth = 1;
        if (isDown){
            ctx.beginPath();
            pts_rect(ctx,e0,e1);
            // ctx.rect(0,0,100,100);
            ctx.stroke();            
        }else{
            ctx.beginPath();
        }
        // console.log(meta);
        window.pt0 = pt0;
        window.meta=meta;
        // window.meta = JSON.parse(JSON.stringify(meta));
        // console.log("[lastEvent]"+lastEvent);
    }
});
Chart.controllers.scatterWithBox = custom;

