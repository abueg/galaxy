define(["mvc/ui/ui-frames"],function(a){return Backbone.View.extend({initialize:function(b){var c=this;b=b||{},this.frames=new a.View({visible:!1}),this.setElement(this.frames.$el),this.buttonActive=b.collection.add({id:"enable-scratchbook",icon:"fa-th",tooltip:"Enable/Disable Scratchbook",onclick:function(){c.active=!c.active,c.buttonActive.set({toggle:c.active,show_note:c.active,note_cls:c.active&&"fa fa-check"}),!c.active&&c.frames.hide()},onbeforeunload:function(){return c.frames.length()>0?"You opened "+c.frames.length()+" frame(s) which will be lost.":void 0}}),this.buttonLoad=b.collection.add({id:"show-scratchbook",icon:"fa-eye",tooltip:"Show/Hide Scratchbook",show_note:!0,visible:!1,onclick:function(){c.frames.visible?c.frames.hide():c.frames.show()}}),this.frames.on("add remove",function(){this.visible&&0==this.length()&&this.hide(),c.buttonLoad.set({note:this.length(),visible:this.length()>0})}).on("show hide ",function(){c.buttonLoad.set({toggle:this.visible,icon:this.visible&&"fa-eye"||"fa-eye-slash"})})},addDataset:function(a){var b=this;require(["mvc/dataset/data"],function(c){var d=new c.Dataset({id:a});$.when(d.fetch()).then(function(){var a={title:d.get("name"),menu:[{icon:"fa fa-chevron-circle-left",tooltip:"Previous in History",onclick:function(){alert("previous")}},{icon:"fa fa-chevron-circle-right",tooltip:"Next in History",onclick:function(){alert("next")}}]},e=_.find(["tabular","interval"],function(a){return-1!==d.get("data_type").indexOf(a)});if(e){var f=new c.TabularDataset(d.toJSON());_.extend(a,{content:function(a){c.createTabularDatasetChunkedView({model:f,parent_elt:a,embedded:!0,height:"100%"})}})}else _.extend(a,{url:Galaxy.root+"datasets/"+d.id+"/display/?preview=True"});b.add(a)})})},addTrackster:function(a){var b=this;require(["viz/visualization","viz/trackster"],function(c,d){var e=new c.Visualization({id:a});$.when(e.fetch()).then(function(){var a=new d.TracksterUI(Galaxy.root),c={title:e.get("name"),type:"other",content:function(b){var c={container:b,name:e.get("title"),id:e.id,dbkey:e.get("dbkey"),stand_alone:!1},d=e.get("latest_revision"),f=d.config.view.drawables;_.each(f,function(a){a.dataset={hda_ldda:a.hda_ldda,id:a.dataset_id}}),view=a.create_visualization(c,d.config.viewport,d.config.view.drawables,d.config.bookmarks,!1)}};b.add(c)})})},add:function(a){if("_blank"==a.target)window.open(a.url);else if("_top"==a.target||"_parent"==a.target||"_self"==a.target)window.location=a.url;else if(this.active)this.frames.add(a);else{var b=$(window.parent.document).find("#galaxy_main");if("galaxy_main"==a.target||"center"==a.target)if(0===b.length){var c=a.url;c+=-1==c.indexOf("?")?"?":"&",c+="use_panels=True",window.location=c}else b.attr("src",a.url);else window.location=a.url}}})});
//# sourceMappingURL=../../maps/layout/scratchbook.js.map