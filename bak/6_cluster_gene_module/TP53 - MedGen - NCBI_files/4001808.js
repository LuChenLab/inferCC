if (typeof (jQuery) != 'undefined') {
    (function ($) {
        $(function () {            
            var min = Math.ceil(1);
            var max = Math.floor(100000);
            var randomNum = Math.floor(Math.random() * (max - min)) + min;
            var surveyUrl = "/projects/Gene/portal/surveys/seqdbui-survey.js?rando=" + randomNum.toString();
            $.getScript(surveyUrl, function () {
                try {
                    ncbi.seqDbUISurvey.init();    
                } catch (err) {
                    console.info(err);
                }
                
            }).fail(function (jqxhr, settings, exception) {
                console.info('Cannot load survey script', jqxhr);
            });;
        });
    })(jQuery);
};
;
jQuery(function(){
    jQuery('#rss_createfeed').bind('click',createRssFeed);
    function createRssFeed (e){
        e.preventDefault();
        var oThis = jQuery(this);
	   	var args = {
            'QueryKey': oThis.data('qk'),
            'Db': oThis.data('db'),
            'RssFeedName': jQuery('#rss_name').val(),
            'RssFeedLimit': jQuery('#rss_results').val(),
            'HID': oThis.data('hid')
        };
        Portal.$send('CreateRssFeed',args);
    }  
});

;
(function($){

    $(function() {    

        var theSearchInput = $("#term");
        var originalTerm = $.trim(theSearchInput.val());
        var theForm = jQuery("form").has(theSearchInput);
        var dbNode = theForm.find("#database");
        var currDb = dbNode.val();
        var sbConfig = {};
        try{
            sbConfig = eval("({" + theSearchInput.data("sbconfig") + "})");
        }catch(e){}
        var defaultSubmit =  sbConfig.ds == "yes";
        var searched = false;
        var dbChanged = null; //since db.change is triggered as a work around for JSL-2067 
        var searchModified = false; //this is used to allow searching when something esle changed on the page with out the term changing
    
        if(!$.ncbi)
            $.extend($,{ncbi:{}});
        if(!$.ncbi.searchbar)
            $.extend($.ncbi,{searchbar:{}});
            
        $.extend($.ncbi.searchbar,
            (function(){
                //*****************private ******************/
               function doSearchPing() {
                   try{
                    var cVals = ncbi.sg.getInstance()._cachedVals;
                    var searchDetails = {}
                    searchDetails["jsEvent"] = "search";
                    var app = cVals["ncbi_app"];
                    var db = cVals["ncbi_db"];
                    var pd = cVals["ncbi_pdid"];
                    var pc = cVals["ncbi_pcid"];
                    var sel = dbNode[0];
                    var searchDB = sel.options[sel.selectedIndex].value;
                    var searchText = theSearchInput[0].value;
                    if( app ){ searchDetails["ncbi_app"] = app.value; }
                    if( db ){ searchDetails["ncbi_db"] = db.value; }
                    if( pd ){ searchDetails["ncbi_pdid"] = pd.value; }
                    if( pc ){ searchDetails["ncbi_pcid"] = pc.value; }
                    if( searchDB ){ searchDetails["searchdb"] = searchDB;}
                    if( searchText ){ searchDetails["searchtext"] = searchText;}
                    ncbi.sg.ping( searchDetails );
                   }catch(e){
                       console.log(e);
                   }
                }
                function getSearchUrl(term){
                    var url = "";
                    if (typeof(NCBISearchBar_customSearchUrl) == "function") 
                            url = NCBISearchBar_customSearchUrl();
                    if (!url) {
                        var searchURI = dbNode.find("option:selected").data("search_uri");
                        url = searchURI ?  searchURI.replace('$',term) : 
                             "/" + dbNode.val() + "/" + ( term !="" ? "?term=" + term : "");
                        }
                    return url;
                }
            
                return {
                    //*****************exposed attributes and functions ******************/
                    'theSearchInput':theSearchInput,
                    'theForm':theForm,
                    'dbNode':dbNode,
                    'searched':searched,
                    'setSearchModified':function() { searchModified = true; },
                    'setSearchUnmodified':function() { searchModified = false; },
                    'searchModified':function(){return searchModified;},
                    'doSearch':function(e){
                           e.stopPropagation();
                           e.preventDefault();
                           //checking for the searched flag is necessary because the autocompelete control fires on enter key, the form submit also fires on enter key
                           if(searched == false){
                               searched = true;
                               theForm.find('input[type="hidden"][name^="p$"]').attr('disabled', 'disabled');
                               //$("input[name]").not(jQuery(".search_form *")).attr('disabled', 'disabled');
                               if (defaultSubmit)
                                   $.ncbi.searchbar.doSearchPing();
                               else {
                                   var term = $.trim(theSearchInput.val());
                                   if (dbChanged || searchModified || term !== originalTerm){
                                       $.ncbi.searchbar.doSearchPing();
                                       var searchUrl = $.ncbi.searchbar.getSearchUrl(encodeURIComponent(term).replace(/%20/g,'+'));
                                       var doPost = (term.length  > 2000) ? true : false; 
                                       if (doPost){
                                           if (e.data.usepjs){
                                               Portal.$send('PostFrom',{"theForm":theForm,"term":term,"targetUrl":searchUrl.replace(/\?.*/,'')});
                                           }
                                           else{
                                               theForm.attr('action',searchUrl.replace(/\?.*/,''));
                                               theForm.attr('method','post');
                                           }
                                       }
                                       else {
                                           window.location = searchUrl;
                                       }
                                   }
                                   else{ //if (term !== originalTerm){
                                       searched = false;
                                   }
                               }
                           }
                    },
                    'onDbChange':function(e){
                         if (dbChanged === null)
                             dbChanged = false;
                         else
                             dbChanged = true;
                         var optionSel = $(e.target).find("option:selected");
                         var dict = optionSel.data("ac_dict");
                         if (dict){
                             //theSearchInput.ncbiautocomplete("option","isEnabled",true).ncbiautocomplete("option","dictionary",dict);
                             theSearchInput.ncbiautocomplete().ncbiautocomplete({
                                    isEnabled: true,
                                    dictionary: dict
                                });
                             theSearchInput.attr("title","Search " + optionSel.text() + ". Use up and down arrows to choose an item from the autocomplete.");
                         }
                         else{
                           theSearchInput.ncbiautocomplete().ncbiautocomplete("turnOff",true);
                           theSearchInput.attr("title", "Search " + optionSel.text());
                         }
                         if (defaultSubmit)
                            theForm.attr('action','/' + dbNode.val() + '/');  
                    },
                    'doSearchPing':function(){
                        doSearchPing();
                    },
                    'getSearchUrl':function(term){
                        return getSearchUrl(term);
                    }
                    
                };//end of return 
             })() //end of the self executing anon
        );//end of $.extend($.ncbi.searchbar
    
         function initSearchBar(usepjs){
            //enable the controls for the back button
            theForm.find('input[type="hidden"][name^="p$"]').removeAttr('disabled');
             if (usepjs)
                 portalSearchBar();
         }
         
        
    
        function portalSearchBar(){
            
            Portal.Portlet.NcbiSearchBar = Portal.Portlet.extend ({
                init:function(path,name,notifier){
                    this.base (path, name, notifier);
                },
                send:{
                    "Cmd":null,
                    "Term":null
                },
                "listen":{
                    "PostFrom":function(sMessage,oData,sSrc){
                        this.postForm(oData.theForm,oData.term,oData.targetUrl);
                    }
                },
                "postForm":function(theForm,term,targetUrl){
                       //console.log('targetUrl = ' + targetUrl);
                       theForm.attr('action',targetUrl);
                       theForm.attr('method','post');
                       this.send.Cmd({
                            'cmd' : 'Go'
                        });
                           this.send.Term({
                            'term' : term
                        });
                        Portal.requestSubmit();
                },
                'getPortletPath':function(){
                    return this.realpath + '.Entrez_SearchBar';
                }
            });
    
        }//portalSearchBar
        


         //portal javascript is required to make a POST when the rest of the app uses portal forms 
         var usepjs = sbConfig.pjs == "yes"; 
         //console.log('sbConfig',sbConfig);
         initSearchBar(usepjs);
         
         dbNode.on("change",$.ncbi.searchbar.onDbChange);
        
        theForm.on("submit",{'usepjs':usepjs},$.ncbi.searchbar.doSearch);
        theSearchInput.on("ncbiautocompleteenter ncbiautocompleteoptionclick", function(){theForm.submit();});
        //a work around for JSL-2067
        dbNode.trigger("change");
        //iOS 8.02 changed behavior on autofocus, should probably check other mobile devices too
        if (sbConfig.afs == "yes" && !/(iPad|iPhone|iPod)/g.test(navigator.userAgent) ){ 
            window.setTimeout(function(){
                try{
                	var x = window.scrollX, y = window.scrollY; // EZ-8676
                	
                    var size= originalTerm.length;
                    if (size == 0 || /\s$/.test(originalTerm))
                        theSearchInput.focus()[0].setSelectionRange(size, size);
                    else
                        theSearchInput.focus().val(originalTerm + " ")[0].setSelectionRange(size+1, size+1);
                        
                    window.scrollTo(x, y);
                }
                catch(e){} //setSelectionRange not defined in IE8
            },1);
        }
        
        //set the query changed flag true after a few seconds, still prevents scripted clicking or stuck enter key
        window.setTimeout(function(){$.ncbi.searchbar.setSearchModified();},2000);
         
     });//End of DOM Ready

})(jQuery);

/*
a call back for the 'Turn off' link at the bottom of the auto complete list
*/
function NcbiSearchBarAutoComplCtrl(){
    jQuery("#term").ncbiautocomplete("turnOff",true);
    if (typeof(NcbiSearchBarSaveAutoCompState) == 'function')
        NcbiSearchBarSaveAutoCompState();
 }

 



;
jQuery(function () {
    Portal.Portlet.Entrez_SearchBar = Portal.Portlet.NcbiSearchBar.extend ({
        init:function(path,name,notifier){
            this.base (path, name, notifier);
            var oThis = this;
            jQuery("#database").on("change", function(){
                oThis.send.DbChanged({'db' : this.value});
            });
        },
        send:{
            "Cmd":null,
            "Term":null,
            "DbChanged":null
        },
        'listen':{
            "PostFrom":function(sMessage,oData,sSrc){
        	    this.postForm(oData.theForm,oData.term,oData.targetUrl);
        	    },
            "ChangeAutoCompleteState": function(sMessage, oData, sSrc) {
        	    this.ChangeAutoCompleteState(sMessage, oData, sSrc);
                },
            'CreateRssFeed':function(sMessage,oData,sSrc){
                this.createRssFeed(sMessage,oData,sSrc);
            },
            'AppendTerm': function(sMessage, oData, sSrc) {
    		    this.ProcessAppendTerm(sMessage, oData, sSrc);
    		},
    		// to allow any other portlet to clear term if needed  
    		'ClearSearchBarTerm': function(sMessage, oData, sSrc) {
    			jQuery("#term").val("");
    		},
    		// request current search bar term to be broadcast  
    		'SendSearchBarTerm': function(sMessage, oData, sSrc) {
    			this.send.Term({'term' : jQuery("#term").val()});
    		}
        },
        'createRssFeed':function(sMessage,oData,sSrc){
            
            var site = document.forms[0]['p$st'].value;
    	   	var portletPath = this.getPortletPath();
    	   	
            try{
                var resp = xmlHttpCall(site, portletPath, 'CreateRssFeed', oData, receiveRss, {}, this);
            }
            catch (err){
                alert ('Could not create RSS feed.');
            }
            function receiveRss(responseObject, userArgs) {
        	    try{
            	    //Handle timeouts 
            	    if(responseObject.status == 408){
            	        //display an error indicating a server timeout
            	        alert('RSS feed creation timed out.');
            	    }
            	    
            	    // deserialize the string with the JSON object 
            	    var response = '(' + responseObject.responseText + ')'; 
            	    var JSONobject = eval(response);
            	    // display link to feed
            	    jQuery('#rss_menu').html(JSONobject.Output,true);
            	    //jQuery('#rss_dropdown a.jig-ncbipopper').trigger('click');
            	    jQuery('#rss_dropdown a.jig-ncbipopper').ncbipopper('open');
            	    //document.getElementById('rss_menu').innerHTML = JSONobject.Output;
                }
                catch(e){
                    alert('RSS unavailable.');
                }
            }
                
        },
        'getPortletPath':function(){
            return this.realpath + '.Entrez_SearchBar';
        },
        "ChangeAutoCompleteState": function(sMessage, oData, sSrc){
            var site = document.forms[0]['p$st'].value;
            var resp = xmlHttpCall(site, this.getPortletPath(), "ChangeAutoCompleteState", {"ShowAutoComplete": 'false'}, function(){}, {}, this);
        },
        "ProcessAppendTerm" : function(sMessage, oData, sSrc){
            var theInput = jQuery("#term");
    	    var newTerm = theInput.val();
    	    if (newTerm != '' && oData.op != ''){
    	        newTerm = '(' + newTerm + ') ' + oData.op + ' ';
    	    }
    	    newTerm += oData.term;
    	    theInput.val(newTerm); 
    	    
    	    theInput.focus();
    	}
    }); //end of Portlet.extend
}); //end of jQuery ready

function NcbiSearchBarSaveAutoCompState(){
    Portal.$send('ChangeAutoCompleteState');
}


;
jQuery(function () {
Portal.Portlet.MedGen_SearchBar = Portal.Portlet.Entrez_SearchBar.extend ({
  
  init: function (path, name, notifier) {
    this.base (path, name, notifier);
  },
  
  /* ######### this is a hack. See detailed comment on same function in base */
  "getPortletPath" : function(){
      return (this.realname + ".Entrez_SearchBar");
  }
});
});


;
(function( $ ){ // pass in $ to self exec anon fn

    // on page ready    
    
        $( 'div.portlet' ).each( function() {

            // get the elements we will need
            var $this = $( this );
            var anchor = $this.find( 'a.portlet_shutter' );
            var portBody = $this.find( 'div.portlet_content' );

            // we need an id on the body, make one if it doesn't exist already
            // then set toggles attr on anchor to point to body
            var id = portBody.attr('id') || $.ui.jig._generateId( 'portlet_content' );
            portBody.attr('id', id );
            anchor.attr('toggles', id );

            // initialize jig toggler with proper configs, then remove some classes that interfere with 
            // presentation
            var togglerOpen = anchor.hasClass('shutter_closed')? false : true; 
            anchor.ncbitoggler({
                isIcon: false,
                initOpen: togglerOpen 
            }).
                removeClass('ui-ncbitoggler-no-icon').
                removeClass('ui-widget');

            // get rid of ncbitoggler css props that interfere with portlet styling, this is hack
            // we should change how this works for next jig release
            anchor.css('position', 'absolute').
                css('padding', 0 );

            $this.find( 'div.ui-helper-reset' ).
                removeClass('ui-helper-reset');

            portBody.removeClass('ui-widget').
                css('margin', 0);

            // trigger an event with the id of the node when closed
            anchor.bind( 'ncbitogglerclose', function() {
                anchor.addClass('shutter_closed');
            });

            anchor.bind('ncbitoggleropen', function() {
                anchor.removeClass('shutter_closed');
            });

        });  // end each loop and end on page ready
})( jQuery );
/*
jQuery(document).bind('ncbitogglerclose ncbitoggleropen', function( event ) {
           var $ = jQuery;
           var eventType = event.type;
           var t = $(event.target);
           
          alert('event happened ' + t.attr('id'));
   
           if ( t.hasClass('portlet_shutter') || false ) { // if it's a portlet
               // get the toggle state
               var sectionClosed = (eventType === 'ncbitogglerclosed')? 'true' : 'false';
               alert ('now call xml-http');

            }
        });
*/

Portal.Portlet.NCBIPageSection = Portal.Portlet.extend ({
	init: function (path, name, notifier){
		this.base (path, name, notifier);
		
		this.AddListeners();
	},
    
	"AddListeners": function(){
        var oThis = this;
        
		jQuery(document).bind('ncbitogglerclose ncbitoggleropen', function( event ) {
            var $ = jQuery;
            var eventType = event.type;
            var t = $(event.target);
            
            // proceed only if this is a page section portlet {
            if ( t.hasClass('portlet_shutter')){
                var myid = '';
                if (oThis.getInput("Shutter")){
                    myid = oThis.getInput("Shutter").getAttribute('id');
                }
    
                // if the event was triggered on this portlet instance
                if (t.attr('id') && t.attr('id') == myid){
                    // get the toggle state
                    var sectionClosed = (eventType === 'ncbitogglerclose')? 'true' : 'false';
                    // react to the toggle event
                    oThis.ToggleSection(oThis.getInput("Shutter"), sectionClosed);
                }
            } // if portlet            
        });
	},
	
	"ToggleSection": function(target, sectionClosed){
	   // if remember toggle state, save the selection and log it
	   if (target.getAttribute('remembercollapsed') == 'true'){
	       this.UpdateCollapsedState(target, sectionClosed);
	   }else {
	       this.LogCollapsedState(target, sectionClosed);
	   }
	},
	
	"UpdateCollapsedState": function(target, sectionClosed){
	    var site = document.forms[0]['p$st'].value;
	    var args = { "PageSectionCollapsed": sectionClosed, "PageSectionName": target.getAttribute('pgsec_name')};
	    // Issue asynchronous call to XHR service
        var resp = xmlHttpCall(site, this.getPortletPath(), "UpdateCollapsedState", args, this.receiveCollapse, {}, this);  
	},
	
	"LogCollapsedState": function(target, sectionClosed){
	    var site = document.forms[0]['p$st'].value;
	    // Issue asynchronous call to XHR service
        var resp = xmlHttpCall(site, this.getPortletPath(), "LogCollapsedState", {"PageSectionCollapsed": sectionClosed}, this.receiveCollapse, {}, this);  
	},
	
	'getPortletPath': function(){
        return this.realname;
    }, 
    
    receiveCollapse: function(responseObject, userArgs) {
    }
	
});
		 
;
Portal.Portlet.SensorPageSection = Portal.Portlet.NCBIPageSection.extend ({
	init: function (path, name, notifier){
		this.base (path, name, notifier);
	}
});

(function( $ ){ // pass in $ to self exec anon fn

    // on page ready
    $( function() {
    
        $( 'div.sensor' ).each( function() {

            // get the elements we will need
            var $this = $( this );
            var anchor = $this.find( 'a.portlet_shutter' );
            var portBody = $this.find( 'div.sensor_content' );

            // we need an id on the body, make one if it doesn't exist already
            // then set toggles attr on anchor to point to body
            var id = portBody.attr('id') || $.ui.jig._generateId( 'sensor_content' );
            portBody.attr('id', id );
            anchor.attr('toggles', id );

            // initialize jig toggler with proper configs, then remove some classes that interfere with 
            // presentation
            var togglerOpen = anchor.hasClass('shutter_closed')? false : true; 
            anchor.ncbitoggler({
                isIcon: false,
                initOpen: togglerOpen 
            }).
                removeClass('ui-ncbitoggler-no-icon').
                removeClass('ui-widget');

            // get rid of ncbitoggler css props that interfere with portlet styling, this is hack
            // we should change how this works for next jig release
            anchor.css('position', 'absolute').
                css('padding', 0 );

            $this.find( 'div.ui-helper-reset' ).
                removeClass('ui-helper-reset');

            portBody.removeClass('ui-widget').
                css('margin', 0);

            // trigger an event with the id of the node when closed
            anchor.bind( 'ncbitogglerclose', function() {
                anchor.addClass('shutter_closed');
            });

            anchor.bind('ncbitoggleropen', function() {
                anchor.removeClass('shutter_closed');
            });

        });  // end each loop          
    });// end on page ready
})( jQuery );
;
Portal.Portlet.GeneSensor = Portal.Portlet.SensorPageSection.extend({
    init: function (path, name, notifier) {
        var oThis = this;
		console.info("Created GeneSensor");
		this.base(path, name, notifier);
		
		/*To add linkpos to species links*/
		jQuery(".speciesline .specieslink").each(function(index){
		    var ref = jQuery(this).attr('ref');
		    ref= ref+"&linkpos="+(index+1);
		    jQuery(this).attr('ref', ref);
		});
    },
    
    send: { 
        'Cmd': null, 
        //'DbChanged': null,
        'LinkName': null,
        //'PresentationChange': null,
        'LastQueryKey': null       
    }, 
    
    listen: {       
        'GeneLink<click>':function (e, target, name) {
            //this.send.DbChanged({'db' : 'gene'});
            this.send.Cmd({'cmd' : 'link'});
            this.send.LinkName({'linkname' : 'gene_pubmed_rif'});
            this.send.LastQueryKey({'qk': target.getAttribute('value')});           
            Portal.requestSubmit(); 
        }        
    }
});
;
Portal.Portlet.Entrez_DisplayBar = Portal.Portlet.extend({

	init: function(path, name, notifier) {
		console.info("Created DisplayBar");
		this.base(path, name, notifier);
		
		// for back button compatibility reset values when page loads
		if (this.getInput("Presentation")){
		    this.setValue("Presentation", this.getValue("LastPresentation"));
		    Portal.Portlet.Entrez_DisplayBar.Presentation = this.getValue("LastPresentation");
		}
		if (this.getInput("Format")){
		    this.setValue("Format", this.getValue("LastFormat"));
		    Portal.Portlet.Entrez_DisplayBar.Format = this.getValue("LastFormat");
		}
		if (this.getInput("PageSize")){
		    this.setValue("PageSize", this.getValue("LastPageSize"));
		    Portal.Portlet.Entrez_DisplayBar.PageSize = this.getValue("LastPageSize");
		}
		if (this.getInput("Sort")){
		    this.setValue("Sort", this.getValue("LastSort"));
		    Portal.Portlet.Entrez_DisplayBar.Sort = this.getValue("LastSort");
		}
		this.ResetDisplaySelections();
		this.ResetSendToSelection();
		
    	jQuery( 
            function(){
        
                var animationTime = jQuery("#sendto2").ncbipopper("option","openAnimationTime");
                var currentCnt = 0;
                var expTimer;
        
                function testPosition(){
                    jQuery(window).trigger("ncbipopperdocumentresize");
                    currentCnt+=10;
                    if (currentCnt<animationTime) {
                        expTimer = window.setTimeout(testPosition,10);
                    }
                }
        
                jQuery("#send_to_menu2 input").on("change click", 
                    function(){
                        currentCnt = 0;
                        if(expTimer) window.clearTimeout(expTimer);
                        testPosition();
                    } 
                );
        
            }
        );
		        
	},
	
	
	send: {
		'Cmd': null, 
		'PageSizeChanged': null,
		'ResetSendTo': null,
		'ResetCurrPage': null,
		'AddUserMessage': null
	},
		
	
	listen: {
		
		/* browser events */
			
		"sPresentation<click>": function(e, target, name){
		    this.PresentationClick(e, target, name); 
		},
		
		"sPresentation2<click>": function(e, target, name){
		    this.PresentationClick(e, target, name); 
		},
		
		"sPageSize<click>": function(e, target, name){	
		    this.PageSizeClick(e, target, name);
		},
		
		"sPageSize2<click>": function(e, target, name){	
		    this.PageSizeClick(e, target, name);
		},
		
		"sSort<click>": function(e, target, name){
		    this.SortClick(e, target, name);
		},
		
		"sSort2<click>": function(e, target, name){
		    this.SortClick(e, target, name);
		},
		
		"SetDisplay<click>": function(e, target, name){
			this.DisplayChange(e, target, name); 
		},
		
		"SendTo<click>": function(e, target, name){
			var sendto = target.value;
            var idx = target.getAttribute('sid') > 10? "2" : "";
			this.SendToClick(sendto, idx, e, target, name); 
		},
		
		"SendToSubmit<click>": function(e, target, name){
		    e.preventDefault();
		    var cmd = target.getAttribute('cmd').toLowerCase();
		    var idx = target.getAttribute('sid') > 10? "2" : "";
			this.SendToSubmitted(cmd, idx, e, target, name); 
		},
		
		/* messages from message bus*/
		
		'ResetSendTo' : function(sMessage, oData, sSrc) {
		    this.ResetSendToSelection();
		}
	
	}, // end listen
	
	
	
	/* functions */
	
	'PresentationClick': function(e, target, name){
		Portal.Portlet.Entrez_DisplayBar.Presentation = target.value;
		Portal.Portlet.Entrez_DisplayBar.Format = target.getAttribute('format');
		this.DisplayChange();
	},
	
	'PageSizeClick': function(e, target, name){ 
		Portal.Portlet.Entrez_DisplayBar.PageSize = target.value;
		this.DisplayChange();
	},
	
	'SortClick': function(e, target, name){
		Portal.Portlet.Entrez_DisplayBar.Sort = target.value;
		this.DisplayChange();
	},
	
	'DisplayChange': function(e, target, name){
	    var submit = false;
	    var extractdb = window.location.pathname.match(/\/([A-Za-z]+)\/?/); 
	    var db = (extractdb[1] && extractdb[1] != '') ? extractdb[1] : "";
	    
	    if (db != '' && getEntrezSelectedItemCount() == 1){
	        //get id, attach db and report, and link	        
	        var URL = '/' + db + '/' + getEntrezSelectedItemList() + '?report=' + Portal.Portlet.Entrez_DisplayBar.Presentation
	        + (Portal.Portlet.Entrez_DisplayBar.Format.toLowerCase() == 'text' ? '&format=text' : '');
	        window.location = URL;
	    }
	    else if (db != '' && getEntrezResultCount() == 1 && window.location.href != ""){   
	        //remove report= from URL and insert new report= into URL
	        if ((window.location.pathname != '' && window.location.pathname.match(/\/[A-Za-z]+\/\w*\d+\w*/))
	            || window.location.href.match(/\/[A-Za-z]+\/??.*term=[^&\s]+/)
	        ){
	            var URL = window.location.href.replace(/&?report=\w+/, "").replace(/\?&/, "?");
	            var hashtagindex = URL.indexOf("#");
	            if (hashtagindex >= 0){
	                URL = URL.substring(0, hashtagindex);
	            }
	            URL += (URL.match(/\?/) ? (URL.match(/\?[^\s]+/) ? "&" : "") : "?") 
	                + "report=" + Portal.Portlet.Entrez_DisplayBar.Presentation
	                + (Portal.Portlet.Entrez_DisplayBar.Format.toLowerCase() == 'text' ? '&format=text' : '');
	            window.location = URL;    
	        }
	        else {
	            submit = true;
	        }
	    }
	    else{
            submit = true;
        }
        
        if (submit){
            this.send.Cmd({'cmd': 'displaychanged'});
            
    	    this.SetPresentationChange(e, target, name);
    	    this.SetPageSizeChange(e, target, name);
    	    this.SetSortChange(e, target, name);
    	    
    	    Portal.requestSubmit();
	    }
	},
	
	'SetPresentationChange': function(e, target, name){
        this.setValue("Presentation", Portal.Portlet.Entrez_DisplayBar.Presentation);
	    this.setValue("Format", Portal.Portlet.Entrez_DisplayBar.Format);
	},
	
	'SetPageSizeChange': function(e, target, name){
	    this.setValue("PageSize", Portal.Portlet.Entrez_DisplayBar.PageSize);
		if (this.getValue("PageSize") != this.getValue("LastPageSize")){
    		//send PageSizeChanged
    		this.send.PageSizeChanged({
    			'size': this.getValue("PageSize"),
                'oldsize': this.getValue("LastPageSize")
    		});	
		}
	},
		
	'SetSortChange': function(e, target, name){
	    if (this.getInput("Sort")){
	        this.setValue("Sort", Portal.Portlet.Entrez_DisplayBar.Sort);
            if (this.getValue("Sort") != this.getValue("LastSort")){
                // ask to reset CurrPage 
    		    this.send.ResetCurrPage();
    		}
    		
    		// set sort in cookie   		
    		var extractdb = window.location.pathname.match(/\/([A-Za-z]+)\/?/); 
    	    var db = (extractdb[1] && extractdb[1] != '') ? extractdb[1] : "";
    	    
    		this.SetSortCookie(Portal.Portlet.Entrez_DisplayBar.Sort, db);
        }    	
	},
		
	'SendToClick': function(sendto, idx, e, target, name) {
		if(sendto.toLowerCase() == 'file'){
			this.SendToFile(sendto, idx);
		}
		else if(sendto.toLowerCase() == 'addtocollections'){
			this.SendToCollections(sendto, idx);
		}
		else if(sendto.toLowerCase() == 'addtoclipboard'){
		    this.SendToClipboard(sendto, idx);
		}
	},
	
	'SendToSubmitted': function(cmd, idx, e, target, name){
	    if (cmd == 'addtobibliography'){
	    	this.SendToBibliographySubmitted(e, cmd, idx, target);
	    }
	    else {
    	    if (cmd == 'file'){
    	         this.SendToFileSubmitted(cmd, idx, target);
    	    }
    	    else if (cmd == 'addtocollections'){
    	    	this.SendToCollectionsSubmitted(cmd, idx, target);
    	    }
    	    this.send.Cmd({'cmd': cmd});
    	    Portal.requestSubmit();
        }       
	},
	
	'ResetSendToSelection': function(){
	    var SendToInputs = this.getInputs("SendTo");
	    for (var j = 0; j < SendToInputs.length; j++){
		    if (SendToInputs[j].checked){
		        SendToInputs[j].checked = false;
			}
		}
	},
	
	'SendToFile': function(name, idx){
	    // generate content
	    var count = this.getItemCount();
		var content = 'Download ' + count + ' items.';
		this.addSendToHintContent(name, idx, content);
	},
	
	'SendToCollections': function(name, idx){
	    // generate content
        var count = this.getItemCount();
        var content= 'Add ';
        var optionNode = document.getElementById("coll_start_option" + idx);
        if (count > Portal.Portlet.Entrez_DisplayBar.CollectionsUpperLimit){
            content += Portal.Portlet.Entrez_DisplayBar.CollectionsUpperLimitText;
            if (optionNode){
            	optionNode.className = '';
            }
        }
        else{
            content += count;
            if (optionNode){
            	optionNode.className = 'hidden';
            }
        }
        content += " items.";
        this.addSendToHintContent(name, idx, content);	
	},
	
	'SendToClipboard': function(name, idx){
	    // generate content
	    var count = this.getItemCount();
        var content= 'Add ';
        if (count > Portal.Portlet.Entrez_DisplayBar.ClipboardLimit){
            content += "the first " + Portal.Portlet.Entrez_DisplayBar.ClipboardLimit;
        }
        else{
            content += count;
        }
        content += " items.";
        this.addSendToHintContent(name, idx, content);
	},
	
	'getItemCount': function(){
	    // ask for selected items count from DbConnector
	    var selectedItemCount = getEntrezSelectedItemCount();
	    if (selectedItemCount > 0){
	        return selectedItemCount;
	    }
	    else{
	        // ask for result count from Entrez_ResultsController
	        return getEntrezResultCount();
	    }
	},
	
	'addSendToHintContent': function(name, idx, content){
	    var hintNode = document.getElementById("submenu_" + name + "_hint" + idx);
	    if (hintNode){
	        hintNode.innerHTML = content;
	        hintNode.className = 'hint';
	    }
	},
	
	'AddSendToSubmitEvent': function(){
	    // add event for SendTo submit button click. 
	    // This call is needed if the position of the submit button node has changed in relation to its parent node. 
        this.addEvent("SendToSubmit", "click", function(e, target, name) {
            var cmd = target.getAttribute('cmd');
            this.SendToSubmitted(cmd, e, target, name); 
        }, false);
    },
    
    'SendToFileSubmitted': function(cmd, idx, target){
         if (this.getInput("FFormat" + idx)){
             this.setValue("FileFormat", this.getValue("FFormat" + idx));
         }
         if (this.getInput("FSort" + idx)){
             this.setValue("FileSort", this.getValue("FSort" + idx));
         }
    },
    
    'SendToCollectionsSubmitted': function(cmd, idx, target){
         if (document.getElementById("coll_start" + idx)){
             document.getElementById("coll_startindex").value = document.getElementById("coll_start" + idx).value;
         }
    },
    
    'SendToBibliographySubmitted': function(e, cmd, idx, target){
        var oThis = this;	
        jQuery.ui.jig.requiresLoginURL = "/account/signin/?inlinelogin=true&popuplogin=true";
        jQuery.ui.jig.requiresLogin( function(name, requiredLogin){ 	 
	        oThis.send.Cmd({'cmd': cmd});
	        oThis.send.AddUserMessage({'type': 'info', 
	                                    'name': 'mybib_processing_msg',
	                                    'msg': 'Adding items to bibliography ...'})
	        Portal.requestSubmit();
	        // Hack to directly make portal submit the page from async function call requiresLogin
	        // By the time asnc function finishes execution portal is not looking for submit request. It's too long after an event firing.
	        // Investigated by Mark Johnson, on 1/17/2019
	        d = Dispatcher.getInstance();
            d.submit();
	    }); 
    },
    
    'ResetDisplaySelections': function(){
        if (this.getInput("Presentation")){
            var selection = this.getValue("Presentation").toLowerCase() + this.getValue("Format").toLowerCase();
            if (document.getElementById(selection)){
                document.getElementById(selection).checked = true;
            }
            // bottom display bar
            if (document.getElementById(selection + "2")){
                document.getElementById(selection + "2").checked = true;
            }
            
        }
        if (this.getInput("PageSize")){
            var selection = 'ps' + this.getValue("PageSize");
            if (document.getElementById(selection)){
                document.getElementById(selection).checked = true;
            }
            // bottom display bar
            if (document.getElementById(selection + "2")){
                document.getElementById(selection + "2").checked = true;
            }
        }
        if (this.getInput("Sort")){
            var selection = this.getValue("Sort") || 'none'; 
            if (document.getElementById(selection)){
                document.getElementById(selection).checked = true;
            }
            // bottom display bar
            if (document.getElementById(selection + "2")){
                document.getElementById(selection + "2").checked = true;
            }
        }
    },
    
    'SetSortCookie': function(sort, db){
	    if (db != ''){
            var d = new Date();
            d.setTime(d.getTime() + (365*24*60*60*1000));
            var expires = "expires="+d.toUTCString();
            
            var newCookie = db + ":" + sort;
            var oldCookie = this.getCookie('entrezSort');
            if (oldCookie != ''){
                if (oldCookie.indexOf(db) != -1){
                    var oldSortVal = oldCookie.substring(oldCookie.indexOf(db));
                    if (oldSortVal.indexOf('&') != -1){
                        oldSortVal = oldSortVal.substring(0, oldSortVal.indexOf('&'));
                    }
                    newCookie = oldCookie.replace(oldSortVal, newCookie);
                }
                else{
                    newCookie = newCookie + "&" + oldCookie;
                }
            } 
            newCookie = "entrezSort=" + newCookie + ";domain=.ncbi.nlm.nih.gov;path=/;" + expires;
            document.cookie = newCookie;
            
		}
    },
    
    // from http://www.w3schools.com/js/js_cookies.asp
    'getCookie': function (cname) {
        var name = cname + "=";
        var ca = document.cookie.split(';');
        console.info("cookie count: " + ca.length);
        for(var i=0; i<ca.length; i++) {
            var c = ca[i];
            while (c.charAt(0)==' ') c = c.substring(1);
            if (c.indexOf(name) == 0) return c.substring(name.length,c.length);
        }
        return "";
    }
	
},
{
    Presentation: '',
    Format: '',
    PageSize: '',
    Sort: '',
    CollectionsUpperLimit: 1000,
	CollectionsUpperLimitText: '1,000',
	ClipboardLimit: 500
});



;
//	BioConcepts_DisplayBar: hemv 05/05/2010
Portal.Portlet.MedGen_DisplayBar = Portal.Portlet.Entrez_DisplayBar.extend({

			init: function(path, name, notifier) {
				this.base(path, name, notifier);
	}
});



;
Portal.Portlet.Entrez_ResultsController = Portal.Portlet.extend({

	init: function(path, name, notifier) {
		console.info("Created Entrez_ResultsController");
		this.base(path, name, notifier);
	},	
		
	send: {
	    'Cmd': null
	},
		
	listen: {
	
	    /* page events */
	    
	    "RemoveFromClipboard<click>": function(e, target, name){
            this.RemoveFromClipboardClick(e, target, name);
	    },
	    
		/* messages */
		
		'Cmd': function(sMessage, oData, sSrc){
		    this.ReceivedCmd(sMessage, oData, sSrc);
		},
		
		'SelectedItemCountChanged' : function(sMessage, oData, sSrc){
		    this.ItemSelectionChangedMsg(sMessage, oData, sSrc);
		},
		
		// currently sent by searchbox pubmed in journals 
		'RunLastQuery' : function(sMessage, oData, sSrc){
			if (this.getInput("RunLastQuery")){
				this.setValue ("RunLastQuery", 'true');
			}
		}
		
	},//listen
	
	'RemoveFromClipboardClick': function(e, target, name){
	    if(confirm("Are you sure you want to delete these items from the Clipboard?")){
	        this.send.Cmd({'cmd': 'deletefromclipboard'});
		    Portal.requestSubmit();  
    	}
	},
	
	// fix to not show remove selected items message when Remove from clipboard was clicked directly on one item
	'ReceivedCmd': function(sMessage, oData, sSrc){
	    if (oData.cmd == 'deletefromclipboard'){
	        Portal.Portlet.Entrez_ResultsController.RemoveOneClip = true;
	    }
	},
	
	'ItemSelectionChangedMsg': function(sMessage, oData, sSrc){
	    // do not show any messages if one item from clipbaord was removed with direct click.
	    if (Portal.Portlet.Entrez_ResultsController.RemoveOneClip){
	        Portal.Portlet.Entrez_ResultsController.RemoveOneClip = false;
	    }
	    else{
    		this.SelectedItemsMsg(oData.count);
    	    this.ClipRemoveMsg(oData.count);
    	}
	},
	
	'SelectedItemsMsg': function(count){
	    SelMsgNode = document.getElementById('result_sel');
	    if (SelMsgNode){
	        if (count > 0){
	            SelMsgNode.className = 'result_sel';
 	            SelMsgNode.innerHTML = "Selected: " + count;
 	        }
 	        else {
 	            SelMsgNode.className = 'none';
 	            SelMsgNode.innerHTML = "";
 	        }
	    }
	},
	
	'ClipRemoveMsg': function(count){
	    ClipRemNode = document.getElementById('rem_clips');
 	    if (ClipRemNode){
 	        if (count > 0){
 	            ClipRemNode.innerHTML = "Remove selected items";
 	        }
 	        else {
 	            ClipRemNode.innerHTML = "Remove all items";
 	        }
 	    }
	},
	
	'ResultCount': function(){
	    var totalCount = parseInt(this.getValue("ResultCount"));
	    totalCount = totalCount > 0 ? totalCount : 0;
	    return totalCount;
	}

},
{
    RemoveOneClip: false
});

function getEntrezResultCount() {
    var totalCount = document.getElementById("resultcount") ? parseInt(document.getElementById("resultcount").value) : 0;
	totalCount = totalCount > 0 ? totalCount : 0;
	return totalCount;
}

;
Portal.Portlet.MedGen_ResultsController = 
	Portal.Portlet.Entrez_ResultsController.extend(
	{
		init: function(path, name, notifier) 
		{
			console.info("Created MedGen_ResultsController");
			this.base(path, name, notifier);
	  }
	});

function MedGenParse(s)
{
  var s1 = ", " + s;
	var ss = s1.replace(/, *(\w+):\s*/g,', "$1": ').replace(/^, */,"").replace(/'/g,'"');
	var obj = jQuery.parseJSON("{" + ss + "}");
	return obj;
};
;
Portal.Portlet.Entrez_Pager = Portal.Portlet.extend ({

    init: function (path, name, notifier) {		
		this.base (path, name, notifier);
    },
   
   
    send: {
        'Cmd': null
    },
   
   
    listen: {
		// page events
		"Page<click>" : function(e, target, name){
			this.send.Cmd({'cmd': 'PageChanged'});
			this.setValue("CurrPage", target.getAttribute('page'));
			Portal.requestSubmit();
		},
		
		"cPage<keypress>" : function(e, target, name){
		    // get page event
		    var event = e || utils.fixEvent (window.event);
		    // if page event was trying to submit 
		    if ((event.keyCode || event.which) == 13) {
		        this.NewPage(event, target);
		    } 			
		},
		
		// messages

		// when pagesize is changed, pager adjusts page number to keep displaying the start 
		// the start item of the initial page
		'PageSizeChanged' : function(sMessage, oData, sSrc) {
			if (this.getInput("CurrPage")){
				var start = (oData.oldsize * (this.getValue("CurrPage") - 1)) + 1;
				var newPage = parseInt((start - 1)/oData.size) + 1;
				this.setValue("CurrPage", newPage);
			}
		},
		
		'ResetCurrPage' : function(sMessage, oData, sSrc) {
			if (this.getInput("CurrPage")){
				this.setValue("CurrPage", '1');
			}
		}

    },
    
    'NewPage': function (event, target){
        // stop event propagation
        event.returnValue = false;
        if (event.stopPropagation != undefined)
            event.stopPropagation ();   
        if (event.preventDefault != undefined)
            event.preventDefault ();    
        
        // get new page info
        newPage = target.value.replace(/,/, ''); // remove comma in page number;
        var npage = parseInt(newPage); 
        var lpage = parseInt(target.getAttribute('last'));
        var cpage = this.getValue("CurrPage");
        
        // check validity of new page
        if (!isNaN(lpage) && newPage != cpage){
             // if the page number entered is not a number or it is a negative number
            if (isNaN(npage) || npage <= 0) { 
                alert("This is not a valid page number: " + newPage); 
                target.value = cpage; 
            } 
            // if the entered value was changed during conversion due to forgiving extras
            else if (npage.toString() != newPage) { 
                alert("This is not a valid page number: " + newPage); 
                target.value = cpage; 
            } 
            // if the entered page is larger than the last page
            else if (npage > lpage) {
                alert("This number is outside the page range: " + newPage); 
                target.value = cpage; 
            }
            else {
                // update the page if entered page number was valid
                this.send.Cmd({'cmd': 'PageChanged'});
                this.setValue("CurrPage", newPage);
                Portal.requestSubmit();  
            } 
        }// end if max is not invalid
           
        return false;
    }
    
});

;
Portal.Portlet.Entrez_LimitsTab = Portal.Portlet.extend ({
  
	init: function (path, name, notifier) 
	{ 
		this.base (path, name, notifier);
		
		var limApplyBtn = jQuery("button[name$=ApplyLimits]");
		if(limApplyBtn[0]){
    		jQuery(document.forms[0]).off("submit").on("submit",function(e){
    		    e.preventDefault();
    		    limApplyBtn.trigger("click")
    		});
		}
	},

    /* If you have to OVERRIDE the send and listen sections, please make sure to include the code in 
    those sections here before your code. Same for init if you need to override it. All other functions
    will be carried over by Portal.Portlet.My_LimitsTab = Portal.Portlet.Entrez_LimitsTab.extend .
    */
    
	send: {
		"Cmd": null,
		"SendSearchBarTerm": null
	},
	
	listen: {
	
	    /* actions from the Limits activated message area */
		
		"ChangeLimits<click>": function(e, target, name) {
		    e.preventDefault();
	        e.stopPropagation();
		    this.ProcessChange (e, target, name);
		},
		
		"RemoveLimits<click>": function(e, target, name) {
		    e.preventDefault();
	        e.stopPropagation();
		    this.ProcessRemove (e, target, name);
		},
		
		"Term": function(sMessage, oData, sSrc) {
		    this.ProcessTerm (sMessage, oData, sSrc);
		},
		
		/* actions from the LimitsPage */
		
		"ClearAllLimits<click>": function(e, target, name){
		    e.preventDefault();
	        e.stopPropagation();
		    this.ProcessClearAll(e, target, name);
		},
		
		"ApplyLimits<click>": function(e, target, name){
		    e.preventDefault();
	        e.stopPropagation();
		    this.ProcessApply(e, target, name);
		},
		
		"CancelLimits<click>": function(e, target, name){
		    e.preventDefault();
	        e.stopPropagation();
		    this.ProcessCancel(e, target, name);
		}
		
	}, //end listen
	
	/* actions from the Limits activated message area */
	
	"ProcessChange" : function(e, target, name){
		Portal.Portlet.Entrez_LimitsTab.WaitingForSearchTerm = true;	
		Portal.Portlet.Entrez_LimitsTab.Link = target.href;
		this.send.SendSearchBarTerm();
	},
	
	"ProcessRemove" : function(e, target, name){
		this.send.Cmd({'cmd': 'removelimits'});
		Portal.requestSubmit();
	},
	
	"ProcessTerm": function(sMessage, oData, sSrc) {
        if (Portal.Portlet.Entrez_LimitsTab.WaitingForSearchTerm) { // make sure ProcessChange was clicked
            Portal.Portlet.Entrez_LimitsTab.WaitingForSearchTerm = false; 

		    if (oData.term){
		        Portal.Portlet.Entrez_LimitsTab.Link += "?term=" + escape(oData.term);
	        }
		    window.location = Portal.Portlet.Entrez_LimitsTab.Link;        
        }       
    },
   
    /* actions from the LimitsPage */
    
    "ProcessCancel" : function(e, target, name){
		history.back();
    },
    
	"ProcessApply" : function(e, target, name){
	    this.send.SendSearchBarTerm();
	    this.send.Cmd({'cmd': 'search'});
		Portal.requestSubmit();
	},
	
	// Implementation will have to change depending on database. All Limit options will have to be cleared.
	"ProcessClearAll" : function(e, target, name){
		this.ClearTagTerms();
	},
	
	"ClearTagTerms": function(){
	    this.getInput("LimitsField").options[0].selected = true;
	},
	
	// Provide implementation for this function if you need to collect some data before the form is submitted.
	// This is necessary if you cannot directly depend on getting form element values from this
	// page into portlet attributes in your Limits portlet after submit
	"beforesubmit": function(){
	    return false;
	}

},

{
    'WaitingForSearchTerm': false,
    'Link': ''
});


// Clear all checkboxes inside target node
function setAll(nodeName, value) {
   if (!document.getElementById) return false;
   var node= document.getElementById(nodeName);

   if (node) {
      var cbs = node.getElementsByTagName("INPUT");
      for (var i = 0; i < cbs.length; i++) {
         var cb = cbs[i];
         if (cb.getAttribute("TYPE").toUpperCase() == "CHECKBOX") {
            cb.checked = value;
         } else {
             cb.value = ""; 
		 }
      }
   }
   return false;
}

;
Portal.Portlet.MedGen_LimitsTab = Portal.Portlet.Entrez_LimitsTab.extend({
	init: function (path, name, notifier) {
	
	
	
		var oThis = this;
		this.base(path, name, notifier);
		
		// update each toggler status
		this.UpdateWholeSelectionInfo('#limit-semantic');
		this.UpdateWholeSelectionInfo('#limit-vocab');
    this.jqDefaultChecked = jQuery('#GTR,#disease_focus');
    function UpdateAllHidden()
    {
      oThis.UpdateAllHidden();
      
    }
		
		// bind event for all checkboxes
		jQuery('div.limboxcontent input:[type=checkbox]').bind('click', function (e) {
			var thisChk = this;
			oThis.UpdateSelectionInfo(thisChk);
		});
		jQuery('.reset-limit').bind('click',
			function(e)
			{
			  var oDiv = jQuery(this).parents('.limboxtitle').parent('.limbox');
			  oThis.MedGenClearSome(oDiv);
			  return false;
			});
    jQuery('button[name*="ClearAllLimits"]').each(
      function(i)
      {
        jQuery(this).attr('title','Set all default limits');
      });
	  
		oThis.jqDefaultChecked.bind('change',
		  function()
		  {
		    oThis.UpdateHidden(this);
		  });
		oThis.UpdateAllHidden();
	},
	
  "UpdateHidden": function(obj)
  {
    var sSelect = obj.checked ? "" : obj.value;
    var sID = '#NOT' + jQuery(obj).attr('id');
    jQuery(sID).each( 
      function(i)
      {
        this.value = sSelect; 
      });
  },
  "UpdateAllHidden": function()
  {
    var oThis = this;
    this.jqDefaultChecked.each(
      function(i)
      {
        oThis.UpdateHidden(this);
      });
  },

  // for handling checkbox click event
	// update the selection status based on the number of checked boxes.
	"UpdateSelectionInfo": function (thisChk) {
		
		thisChk = jQuery(thisChk);		
		var smallBox = thisChk.parents('div.ui-ncbitoggler');
		
		console.info(smallBox);
		
		// grab all selcted checkboxes
		var selChkBoxes = smallBox.find('input:[type=checkbox][checked]');
		
		console.info(selChkBoxes);
		
		// grab the immediate h5, and update the status
		var h5 = jQuery(smallBox.parent('div').prev('h5')[0]);
		
		if (h5) {
			var numSel = h5.find('span.num-sel');
			numSel.text(selChkBoxes.length)
		}
	},
	
	// Given a small subset of vocaublary
	// update the selection status based on the number of checked boxes.
	"UpdateWholeSelectionInfo": function (bigBoxId) {
	
		console.log("Executed UpdateWholeSelectionInfo")
			
		var smallBox = jQuery(bigBoxId + " .limboxcontent div.ui-ncbitoggler");
		
		console.info(smallBox);
		
		smallBox.each(function (k, v) {
			var oThis = jQuery(this);
			// grab all checkboxes that are checked
			// and update status
			var chkBoxes = oThis.find('input:[type = checkbox][checked]');
			var h5 = jQuery(oThis.parent('div').prev('h5')[0]);
			
			if (h5) {
				var numSel = h5.find('span.num-sel');
				numSel.text(chkBoxes.length)
			}
		});
		
		
		// now show all info
		jQuery(bigBoxId + ' div.limboxcontent  span.label-info').fadeIn();
		jQuery(bigBoxId + ' div.limboxcontent  span.num-sel').fadeIn();
		
	},
	
	listen: {
		"btnSemClearAll<click>": function (e, target, name) {
			setAll('limit-semantic', false);
			this.UpdateWholeSelectionInfo('#limit-semantic');
		},
		"btnSemSelAll<click>": function (e, target, name) {
			jQuery('div#limit-semantic input:[type=checkbox]').attr('checked', true);
			this.UpdateWholeSelectionInfo('#limit-semantic');
		},
		
		"btnSelAll<click>": function (e, target, name) {
			jQuery('div#limit-vocab input:[type=checkbox]').attr('checked', true);
			this.UpdateWholeSelectionInfo('#limit-vocab');
		},
		
		"btnSelUMLS<click>": function (e, target, name) {
			
			jQuery('div#limit-vocab input:[type=checkbox]').attr('checked', true);
			jQuery('div#limit-vocab input:[type=checkbox].owned-by-myncbi').attr('checked', false);
			
			this.UpdateWholeSelectionInfo('#limit-vocab');
		},
		
		"btnSelNCBI<click>": function (e, target, name) {
			
			jQuery('div#limit-vocab input:[type=checkbox]').attr('checked', false);
			jQuery('div#limit-vocab input:[type=checkbox].owned-by-myncbi').attr('checked', true);
			
			this.UpdateWholeSelectionInfo('#limit-vocab');
		},
		
		
		"btnClearAll<click>": function (e, target, name) {
			setAll('limit-vocab', false);
			this.UpdateWholeSelectionInfo('#limit-vocab');
		},
		
		/* actions from the Limits activated message area */
		"ChangeLimits<click>": function (e, target, name) {
		    e.preventDefault();
	        e.stopPropagation();
			this.ProcessChange(e, target, name);
		},
		
		"RemoveLimits<click>": function (e, target, name) {
		    e.preventDefault();
	        e.stopPropagation();
			this.ProcessRemove(e, target, name);
		},

//    "DefaultLimits<click>": function (e, target, name) {
//      this.ProcessDefault(e, target, name);
//    },

		"Term": function (sMessage, oData, sSrc) {
			this.ProcessTerm(sMessage, oData, sSrc);
		},
		
		/* actions from the LimitsPage */
		"ClearAllLimits<click>": function (e, target, name) {
		    e.preventDefault();
	        e.stopPropagation();
			this.ProcessClearAll(e, target, name);
		},
    "DefaultAllLimits<click>": function (e, target, name) {
      this.MedGenDefaultAll();
    },
		
		"ApplyLimits<click>": function (e, target, name) {
		    e.preventDefault();
	        e.stopPropagation();
			this.ProcessApply(e, target, name);
		},
		
		"CancelLimits<click>": function (e, target, name) {
		    e.preventDefault();
	        e.stopPropagation();
			this.ProcessCancel(e, target, name);
		}
	},
	// end listen

//  "ProcessDefault" : function(e, target, name){
//    this.send.Cmd({'cmd': 'defaultlimits'});
//    Portal.requestSubmit();
//  },

	"ProcessDefaultAll" : function(e, target,name){
    this.MedGenDefaultAll();
	},
	"ProcessClearAll": function (e, target, name) {
		this.base(e, target, name);
    this.MedGenClearAll();
	},
	
	"BioTermsClearAll": function () {
		setAll('limit-vocab', false);
		setAll('limit-semantic', false);
		return false;
	},
	"MedGenClearAll": function()
	{
		this.MedGenClearSome('body',false);
    this.UpdateAllHidden();
	},
  "MedGenDefaultAll": function()
  {
    this.MedGenClearSome('body',true);
    this.UpdateAllHidden();
  },
	"MedGenClearSome": function (which,bDefault) {
		jQuery('input[type="checkbox"]',which).each(
			function(i)
			{
		        this.checked = bDefault && jQuery(this).hasClass('defaultChecked');
			});
		jQuery('select[name*="LimitsField"]',which).each(
			function(i)
			{
			  	this.selectedIndex = 0;
			});
		return false;
	},
  "MedGenDefaultSome": function(which)	
  {
    this.MedGenClearSome(which,true);
  }
});




;
Portal.Portlet.Entrez_Messages = Portal.Portlet.extend({

	init: function(path, name, notifier) {
		this.base(path, name, notifier);
		
		this.setMsgAreaClassName();
	},
	
	listen: {
	   /* messages from message bus*/
		
		'AddUserMessage' : function(sMessage, oData, sSrc) {
		    // create new message node
		    var msgnode = document.createElement('li');
		    if (oData.type != ''){
		        msgnode.className = oData.type + ' icon'; 
		    }
		    if (oData.name != ''){
		        msgnode.id = oData.name; 
		    }
		    msgnode.innerHTML = "<span class='icon'>" + oData.msg + "</span>";
		    
		    // add new node as first message in message block (not ads that look like messages)
		    var parent = document.getElementById('msgportlet');
		    if (parent){
    		    var oldnode = document.getElementById(oData.name);
    		    if (oldnode){
    		        parent.removeChild(oldnode);
    		    }
    		    var firstchild = parent.firstChild;
    	        if (firstchild){
                    parent.insertBefore(msgnode, firstchild);
                }
                else{
                    parent.appendChild(msgnode);
                }
                this.setMsgAreaClassName('true');
            }
            //if there was no ul, create one, then insert the li
            else {
                var msgarea = document.getElementById('messagearea');
                if (msgarea){
                    var msgportlet = document.createElement('ul');
                    msgportlet.className = 'messages';
                    msgportlet.id = 'msgportlet';
                    msgportlet.appendChild(msgnode);
                    if (msgarea.firstChild){
                         msgarea.insertBefore(msgportlet, msgarea.firstChild);
                    }
                    else{
                        msgarea.appendChild(msgportlet);
                    }
                    this.setMsgAreaClassName('true');
                }
            }
		},
		
		'RemoveUserMessage' : function(sMessage, oData, sSrc) {
		    var msgnode = document.getElementById(oData.name);
		    if (msgnode){
		        var parent = document.getElementById('msgportlet'); 
		        if (parent){
    		        parent.removeChild(msgnode);
    		        this.setMsgAreaClassName();
    		        // if the parent ul has no children then remove the parent
    		        if (parent.firstChild){}
    		        else {
    		            if (document.getElementById('messagearea')) {
    		                document.getElementById('messagearea').removeChild(parent);
    		            }
    		        }
    		    }
		    }
		}
	}, // end listen
	
	'setMsgAreaClassName' : function(hasMsg){
        var msgarea = document.getElementById('messagearea');
	    if (msgarea){
	        var msgclass = "empty";
	        
    	    // if a message was added, hasMsg is set to true at call time to avoid checks. 
    	    // by default, hasMsg is false.
    	    if (hasMsg == 'true'){
    	        msgclass = "messagearea";
    	    }
    	    else if (msgarea.getElementsByTagName('li').length > 0){
                msgclass = "messagearea"; 
        	}
        	
            msgarea.className = msgclass;
        }
	} // end setMsgAreaClassName
});
		
		
;
Portal.Portlet.Entrez_RVBasicReport = Portal.Portlet.extend({
	
	init: function(path, name, notifier) {
		console.info("Created report portlet");
		this.base(path, name, notifier);
	},
	
	send: {
		'ItemSelectionChanged': null,
		'ClearIdList': null,
		'Cmd': null
	},
	
	listen: {
		"uid<click>" : function(e, target, name){
		    this.UidClick(e, target, name);
		},
		
		"RemoveClip<click>" : function(e, target, name){
		    this.ClipRemoveClick(e, target, name);              
		}
	},
	
	'UidClick': function(e, target, name){	
		this.send.ItemSelectionChanged( { 'id': target.value,
		                                  'selected': target.checked });
	},
	
	'ClipRemoveClick': function(e, target, name){
	    this.send.ClearIdList();
		this.send.Cmd({'cmd': 'deletefromclipboard'});
		this.send.ItemSelectionChanged( { 'id': target.getAttribute('uid'),
		                                  'selected': true });
		Portal.requestSubmit();
	}
});
   

;
Portal.Portlet.MedGen_RVDocSum = Portal.Portlet.Entrez_RVBasicReport.extend({
	init: function(path, name, notifier) {
		this.base(path, name, notifier);
	}
});

//jQuery('body').ready(
// function()
// {   
//	jQuery('input','.rprtnum').css('display','none');
// });
 	

;
//
//  MoreLess - jQuery component to break up either text or two child nodes
//    into a 'brief' or 'verbose' display.
//
//  METHOD 1:  Split By Child Nodes
//
//  For splitting child nodes, 
//    the parameter 'nodeBefore' specifies the number of nodes
//    to appear BEFORE the split.  
//    After the split, there are two nodes,
//      the first node contains a link which allows expansion,
//      the second node contains the remainder of child nodes after
//      the count specified by 'nodeBefore'
//      and each of these two nodes has a link 
//      appended which hides its parent and shows its sibling.  
//      For example:
//
//  <span class="jig-moreless" jig-options="nodeBefore: 1">
//    <span>Main text<span>
//    <span>Alternate text</span>
//  </span>
//
//  After initialization it would be 
//  <span class="jig-moreless" data-jigconfig="nodeBefore: 1">
//    <span style="display:inline;">Main text<a> More &gt;</a><span>
//    <span style="display:none;">Alternate text<a> &lt; </a></span>
//  </span>
//
//  Each of the links above will hide its parent <span> and show its parent's
//  sibling
//
//
//  METHOD 2: inner HTML is text only
//
//  For a node whose child is a long string of text, that text will be split
//  and <span> nodes will be created to provide the functionality above.  For example.
//
//  <span class="jig-moreless">This is a long line of text.  It it also very boring and
//    for that reason you shouldn't show all of it.</span>
//
//  After Initialiazation:
//
//  <span class="jig-moreless">This is a long line 
//      <span style="display: inline;">...<a> More &gt;</a></span>
//      <span style="display: none;">of text.  It it also very boring and
//    for that reason you shouldn't show all of it.<a> &lt; </a></span>
//
//   note that the default length of the abbreviated text is 20 characters,
//
//   options for the attribute jig-config are as follows:
//     nodeBefore: number of child nodes before the split.  
//          Child nodes can be text or elements.  This should NOT
//          be used with the 'shortLength' option described below
//     moreText: 'text string inside of the first hyperlink'
//          default is 'More >'
//     lessText: 'text string inside of the second hyperlink'
//          default is '  <'
//     moreTitle: 'title or pop-up text for the first hyperlink'
//          default is 'Show more information'
//     lessTitle: 'title or pop-up text for the last second hyperlink'
//          default is 'Show less information'
//     shortLength: length of the text to be displayed when not all is selected
//          (integer) default is 20.  This should NOT
//          be used with the 'nodeBefore' option described below
//     ellipses: text to show before the 'moreText' link.  The default is 
//          empty if using the nodeBefore parameter and '...' otherwise
//     class: contents of the class attribute for the new <a> elements
//
jQuery(document).ready(
function()
{
  // constants
  
  var WIDGET_CLASS = 'jig-moreless';
  var DEFAULT_SHORT_LENGTH = 20;
  var MIN_SHORT_LENGTH = 10;
  var CSPACE = String.fromCharCode(160);
  var MORE_TEXT =  CSPACE + "More" + CSPACE + ">";
  var LESS_TEXT =  CSPACE + "<";
  var MORE_TITLE = 'Show more information';
  var LESS_TITLE = 'Show less information';
  var SELECTOR = '>*:lt(2)';

  // utility functions
  /*
  function HtmlEncode(s)
  {
    var sRtn = String(s);
    sRtn = sRtn.replace(/\&/g,'&amp;');
    sRtn = sRtn.replace(/\</g,'&lt;');
    sRtn = sRtn.replace(/\>/g,'&gt;');
    sRtn = sRtn.replace(/\"/g,'&quot;');
    return sRtn;
  } */
  function BeginsLower(x)
  {
    var rtn = (typeof(x) == 'string') &&
      (x.length > 0);
    if(rtn)
    {
      var c = x.charAt(0);
      rtn = (c >= 'a') && (c <= 'z');
    }
    return rtn;
  };

  // Setup()

  function Setup($)
  {
    $.fn.MoreLess = function(arg1,arg2,arg3)
    {
      var MoreLess = 
      {
        FunctionSetup: function(obj)
        {
          var sName;
          for (sName in this)
          {
            if(BeginsLower(sName))
            {
              var thisAttr = this[sName];
              var b = (typeof(thisAttr) === 'function');
              b && (obj[sName] = thisAttr);
            }
          }
        },
        LinkShow: function(b)
        {
          jq1 = jQuery(this).parent();
          jqML = jq1 ? jq1.parent() : null;
          jqML && jqML.MoreLess().show(b);
          return false;
        },
        LinkShowLess: function()
        {
          MoreLess.LinkShow.apply(this,[0]);
          return false;
        },
        LinkShowMore: function() 
        {
          MoreLess.LinkShow.apply(this,[1]);
          return false;
        },
        options: function(a,b)
        {
          var sType = typeof(a);
          var rtn = this;
          if(sType == 'object')
          {
            var obj = this.data('options');
            if(!obj)
            {
              this.data('options',a);
            }
            else
            {
              jQuery.extend(obj,a);
            }
          }
          else if((sType == 'string') || (sType == 'number'))
          {
            var sTypeb = typeof(b);
            if((sTypeb == 'string') || (sTypeb == 'number'))
            {
              var obj = {};
              obj[a] = b;
              rtn = MoreLess.options.apply( this,[obj] );
            }
            else
            {
              var obj = this.data('options');
              rtn = obj && obj[a];
            }
          }
          return rtn;
        }, // options:
        show: function(b)
        {
          var ndx = b ? 0 : 1;
          var jqLess = jQuery(this).data('less');
          var jqMore = jQuery(this).data('more');
          if(jqLess && jqMore)
          {
            var a = ['none','inline'];
            jqLess.css('display',a[ndx]);
            jqMore.css('display',a[1 - ndx]);
          }
        }, // show
        showMore: function()
        {
          this.show(1);
        }, // showMore
        showLess: function()
        {
          this.show(0);
        },  // showLess
        init: function(options)
        {
          var jq = jQuery(this);
          var jqChildren = jq.children();
          var jqContents = jq.contents();
          var nLen = jqChildren.length;
          var nContentLength = jqContents.length;
          var nCutLen = DEFAULT_SHORT_LENGTH;
          var nNodeBefore = 0;
          var tagSpan1 = false;
          var tagSpan2 = false;
          var jqSpan1 = false;
          var jqSpan2 = false;
          var sEllipses = false;

          function SetupSpan()
          {
            tagSpan1 = document.createElement('span');
            tagSpan2 = document.createElement('span');
            jqSpan1 = jQuery(tagSpan1);
            jqSpan2 = jQuery(tagSpan2);
            (typeof(sEllipses) == 'string') &&
              sEllipses.length &&
              jqSpan1.text(sEllipses);
          }
          function AddSpan2(jqa,n)
          {
            var nLen = jqa.length;
            for(i = n; i < nLen; i++)
            {
              jqSpan2.append(jQuery(jqa[i]).detach());
            }
          }


          var sText = jq.text();
          var nTextLen = sText.length;
          (typeof(options) == 'object') || (options = {});
          var sOptions = jq.attr('data-jigconfig') || '';
          if(sOptions.length)
          {
            var obj = MedGenParse(sOptions);
            obj && jQuery.extend(options,obj);
          }
          jq.data('options',options);

          nNodeBefore = Number(options.nodeBefore);
          nCutLen = options.shortLength || DEFAULT_SHORT_LENGTH;
            (isNaN(nCutLen) || (nCutLen < MIN_SHORT_LENGTH)) &&
                      (nCutLen = DEFAULT_SHORT_LENGTH);
          sEllipses = options.ellipses;


          if(nNodeBefore >= 0)
          {
              if(nContentLength > nNodeBefore)
              {              
                SetupSpan();
              }
              jqContents.each(function(i)
              {
                if(this.nodeType && (this.nodeType == 1))
                {
                    var jqThis = jQuery(this);
                    (jqThis.css('display') == 'none') &&
                      jqThis.css('display','inherit');
                }
              });
              AddSpan2(jqContents,nNodeBefore);
          }
          else if(nTextLen <= nCutLen)
          {
            jq = false; // not enough text to split
          }
          else
          {
            // loop through contents until a break point is found
            var nCumLength = 0;
            var sTempText;
            var nTempLen;
            var jqc;
            var i;
            var nNodeType;
            (typeof(sEllipses) == 'string') || (sEllipses = '...')
            var nLenEllipses = sEllipses.length;
            for(i = 0; i < nContentLength; i++)
            {
              domc = jqContents[i];
              nNodeType = domc.nodeType;
              sTempText = jQuery(domc).text();
              nTempLen = sTempText.length;
              if((nCumLength + nTempLen) <= (nCutLen - nLenEllipses))
              {
                nCumLength += nTempLen;
              }
              else if(nNodeType == 1)
              {
                // element
                SetupSpan();
                AddSpan2(jqContents,i);
                i = nContentLength;
              }
              else
              {
                // we are in a text node - great
                var nUse = nCutLen - nCumLength - nLenEllipses;
                var sText1;
                if(nUse < 0) 
                {
                  (nUse = 0);
                  sText1 = ''
                }
                else
                {
                  sText1 = sTempText.substr(0,nUse);
                }
                sText2 = sTempText.substr(nUse);
                jQuery(jqContents[i]).remove();
                jq.append(document.createTextNode(sText1));
                SetupSpan();
                jqSpan2.text(sText2);
                AddSpan2(jqContents,i+1);
                i = nContentLength;
              }
            }
          }
          if(jqSpan1)
          {
            var moreText = options.moreText || MORE_TEXT;
            var lessText = options.lessText || LESS_TEXT;
            var moreTitle = options.moreTitle || MORE_TITLE;
            var lessTitle = options.lessTitle || LESS_TITLE;
            var sClass = options['class'] || '';
            jq.data('less',jqSpan1);
            jq.data('more',jqSpan2);
            var elem1 = document.createElement('a');
            var elem2 = document.createElement('a');
            var jElem1 = jQuery(elem1);
            var jElem2 = jQuery(elem2);
            jElem1.attr('href','#').attr('onclick','return false;');
            jElem2.attr('href','#').attr('onclick','return false;');
            lessTitle.length &&
              jElem2.attr('title',lessTitle);
            moreTitle.length &&
              jElem1.attr('title',moreTitle);
            if(sClass.length)
            {
              jElem1.attr('class',sClass);
              jElem2.attr('class',sClass);
            }
            jElem1.bind('click',MoreLess.LinkShowMore);
            jElem2.bind('click',MoreLess.LinkShowLess);
            jElem1.text(moreText);
            jElem2.text(lessText);
            jqSpan1.append(elem1);
            jqSpan2.append(elem2);
            jq.attr('jig-initialized',1);
            var bDisplay1 = (jqSpan1.css('display') != 'none');
            var bDisplay2 = (jqSpan2.css('display') != 'none');
            if(bDisplay1 == bDisplay2)
            {
              MoreLess.show.apply(this,[options.showAll ? 1 : 0]);             
            }
            jq.append(jqSpan1);
            jq.append(jqSpan2);
          }
          else
          {
            // this should never happen
            jq = false;
          }
          return jq;
        }  // init:
      };
      var rtn = false;
      var fnc = false;
      if(this.attr('jig-initialized'))
      {
        rtn = jQuery(this);
      }
      else
      {
        rtn = MoreLess.init.apply(this,[arg1]);
        rtn &&
          rtn.attr('jig-initialized',1);
      }
      if(rtn)
      {
        MoreLess.FunctionSetup(rtn);
        if(BeginsLower(arg1) && 
          (typeof(fnc = MoreLess[arg1]) == 'function')
          )
        {
          var rtnNext = fnc.apply(this,[arg2,arg3]);
          rtn = rtnNext;
        }
        else if(typeof(arg1) == 'object')
        {
          MoreLess.options.apply(this,[arg1]);
        }
      }
      else
      {
        jQuery(this).removeClass(WIDGET_CLASS);
      }
      return rtn;
    };
    
  };
  
  // setup here
  
  Setup(jQuery);
  
  jQuery('.' + WIDGET_CLASS).each(
    function(i)
    {
      jQuery(this).MoreLess();
    });
  
});

;

Portal.Portlet.Entrez_Filters = Portal.Portlet.extend({
	
	init: function(path, name, notifier) {
		console.info("Created FilterTab");
		this.base(path, name, notifier);
	},
	
	send: {
		'Cmd': null,
		'AppendTerm': null,
		'ClearIdList': null,
		'SendSearchBarTerm' : null,
		'SetTimelineFilter':null
	},
	
	listen: {
			//browser events
		"Filter<click>" : function(e, target, name){
		     this.ProcessFilterClick(e, target, name);
		},
		
		"Pin<click>" : function(e, target, name){
             this.ProcessPinClick(e, target, name);
		},

			// messages
		// back button fix
		'Cmd' : function(sMessage, oData, sSrc){
		     this.CheckCmd(sMessage, oData, sSrc);
		}

	},
	
	"ProcessPinClick" : function(e, target, name){
	    // append filter to search bar term
		newTerm = '\"' + target.getAttribute('filter') + '\"[Filter]';
		this.send.AppendTerm({'op' : 'AND', 'term': newTerm});
		/*
		//for back button compatibility, clear any selected ids.
		this.send.ClearIdList();
		// ask search bar to send the current term, after append. Then send search cmd and request page submit.
		this.send.SendSearchBarTerm();
		this.send.Cmd({'cmd': 'search'});
		Portal.requestSubmit();
		*/
	},
	
	"ProcessFilterClick" : function(e, target, name){
	    this.send.Cmd({'cmd': 'FilterChanged'});
	    this.send.SetTimelineFilter({'TimelineYear':''});
		this.setValue("CurrFilter", target.getAttribute('filter'));
		Portal.requestSubmit();
	},
	
	"CheckCmd" : function(sMessage, oData, sSrc){
	    if (oData.cmd != 'FilterChanged'){
			if(this.getValue("CurrFilter") != this.getValue("LastFilter")){
				this.setValue("CurrFilter", this.getValue('LastFilter'));
				console.info("CurrFilter changed to: " + this.getValue('CurrFilter'));
			}
		}
	}
	
});
   


;
Portal.Portlet.MedGen_FiltersPortlet = 
	Portal.Portlet.Entrez_Filters.extend(
	{
		init: function(path, name, notifier) 
		{
			console.info("Created MedGen_ResultsController");
			this.base(path, name, notifier);
	  }
	});


;
Portal.Portlet.RelatedDataLinks = Portal.Portlet.NCBIPageSection.extend({
    init: function (path, name, notifier) {
        var link_descr;
        console.info("Created RelatedDataLinks Ad");
        this.base(path, name, notifier);
        this.initializeControls();
    },
      
    send: { 
        'SendSavedUidList': null
    }, 
    
    listen: {       
        'rdDatabase<change>' : function (e, target, name) {
            this.setSelectButton();
            this.makeXmlHttpCall();
        },
        
        'rdFind<click>':function (e, target, name) {
            e.preventDefault();
	        e.stopPropagation();
            this.SendLink(e, target, name); 
        },
        
        'rdLinkOption<change>':function (e, target, name) {
            this.SetDescription(e, target, name);
        },
        
        //message from DbConnector with selectedIds
        'newUidSelectionList': function(sMessage, oData, sSrc){
            Portal.Portlet.RelatedDataLinks.selectedIdList = oData.list;
        }
    },    
    
    'getPortletPath' : function(){
        return (this.realname + ".NCBIPageSection");
    },   
    
    'responder' : function (responseObject, userArgs) {
        var seldb = document.getElementById('rdDatabase').selectedIndex;
        if(seldb!=0){
            //use "try" so we can gracefully handle errors
            try {
                // Handle timeouts
                if (responseObject.status == 408) {
                    
                    //display an appropriate error message
                }
                
                //convert the string response into a JavaScript Object
                var resp = responseObject.responseText;
                
                console.debug('This is what was returned from the portlet');
                console.info(resp);
                //why don't you take a look at this in firebug?
                
                resp = '(' + resp + ')';
                json_obj = eval(resp);
                
                console.debug('This is the object that we created from the portlet response');
                console.info(json_obj);
                //now look at what it is.
                
                var link_name = json_obj.response;
                link_name = link_name.split(',');
                
                var link_disp_name = json_obj.response_disp;
                link_disp_name = link_disp_name.split(';');
                
                link_descr = json_obj.response_descr;
                link_descr = link_descr.split('||');
                if(link_descr[0]!=''){
                    document.getElementById('rdDescr').innerHTML = link_descr[0];
                    document.getElementById('rdDescr').style.display = "block";
                }
                
                var link_count = link_name.length;
                if(link_count>0){
                    for(var countr=0;countr<link_count;countr++){
                        if(countr==0)
                            document.getElementById('rdLinkOption').options[countr] = new Option(link_disp_name[countr], link_name[countr], true, false);
                        else
                            document.getElementById('rdLinkOption').options[countr] = new Option(link_disp_name[countr], link_name[countr], false, false);
                    }
                }
                
                if(link_count>1)
                    document.getElementById('rdOption').style.display = "block";
                
                jQuery("#rdFind").ncbibutton("enable");
                document.getElementById('rdDescr').style.display = "block";
                
                //Now we will update the display, and change the second input.
               /* this .setValue('output', json_obj.response);*/
                
                
                
                //now catch any errors that may have occured
            }
            catch (e) {
                //display an appropriate error.  Remember, user-friendly messages!
                alert("Please refresh the page and try again. (" + e + ")");
            }
        }
    }, 
    
    'SendLink': function(e, target, name){ 
        window.location = "/" + this.getValue("rdDatabase") + "?linkname=" +  this.getValue("rdLinkOption")
         + (Portal.Portlet.RelatedDataLinks.selectedIdList != '' ? 
            ("&from_uid=" + Portal.Portlet.RelatedDataLinks.selectedIdList) : 
            (jQuery('#rdqk') && jQuery('#rdqk').val() != '' ? "&querykey=" + jQuery('#rdqk').val() : ""));
        
    },
    
    'initializeControls':function(){
        document.getElementById('rdDatabase') .selectedIndex = 0;
        // Resetting Database select on page load/reload
        this.setSelectButton();
    },
    
    'setSelectButton':function(){
        document.getElementById('rdOption').style.display = "none";
        document.getElementById('rdDescr').style.display = "none";
        document.getElementById('rdFind').disabled = true;
        
        this.deleteOption("rdLinkOption");
    },
    
    /*
    FUNCTION NAME: deleteOption
    Delete all the current options from the required drop down menu
    */
    'deleteOption':function(selectbox){
        while(document.forms[0].elements[selectbox].childNodes.length>0) {
            document.forms[0].elements[selectbox].removeChild(document.forms[0].elements[selectbox].childNodes[0])
        }
    },
    
    'SetDescription':function(e, target, name){
        var selOpt = document.getElementById('rdLinkOption').selectedIndex;
        document.getElementById('rdDescr').innerHTML = link_descr[selOpt];
    },
    
    'makeXmlHttpCall':function(){
        var dbto = this.getValue("rdDatabase");
    
            var siteName = document.forms[0][ 'p$st' ].value;
            
            var portletPath = this .realname;
            
            var actionName = 'XMLHTTPhandler';
            
            var args = {
                'related_data_db' : dbto,
                'Db' : document.getElementById('DbName').value
            };
            
            var callback = this .responder;
            
            var userArgs = {
            };
            
            var oThis = this;
            
            try {
                var response = xmlHttpCall(siteName, portletPath, actionName, args, callback, userArgs, oThis);
            }
            catch (err) {
                alert('The following error has occured: ' + err);
            }
            
    }
},
{
	selectedIdList: ''
});
;
Portal.Portlet.Discovery_SearchDetails = Portal.Portlet.NCBIPageSection.extend ({
	init: function (path, name, notifier){
		this.base (path, name, notifier);		
	},
	
	listen: {	    
	    "SearchDetailsTerm<keypress>": function(e, target, name) {
			var event = e || utils.fixEvent (window.event);
			if ((event.keyCode || event.which) == 13) {
			    // Emulate button click
			    this.SearchDetailsTermPress(event, e, target, name);
			}
		},
	    
        "SearchDetailsQuery<click>":  function(e, target, name) {       
		     this.SearchDetailsQueryClick(e, target, name);
		}		
	},
	
	'getPortletPath' : function(){	    
        return (this.realname + ".NCBIPageSection");
    },   
	
	"SearchDetailsTermPress" : function(event,e, target, name){
    	event.returnValue = false;
    	if (event.stopPropagation != undefined)
              event.stopPropagation ();   
    	if (event.preventDefault != undefined)
              event.preventDefault ();
              
    	this.ProcessSearch (target,e);
    	return false;
	},
	
	"SearchDetailsQueryClick": function(e, target, name){
	    this.ProcessSearch (target,e);
	},
	
	"ProcessSearch": function(target,e){
	    e.preventDefault();
	    e.stopPropagation();
	    if (this.getValue('SearchDetailsTerm') != ''){
	        window.location = "/" + this.getInput('SearchDetailsTerm').getAttribute('db') + "?term=" 
    		 + escape(this.getValue('SearchDetailsTerm')) + "&cmd=DetailsSearch";
    	}
    	else{
    	    alert ('There is no term in the Query Translation box to search.');
    	}
	}
	
});
;
(function( $ ){ // pass in $ to self exec anon fn
    // on page ready
    $( function() {
        $('li.ralinkpopper').each( function(){
            var $this = $( this );
            var popper = $this;
            var popnode = $this.find('div.ralinkpop');
            var popid = popnode.attr('id') || $.ui.jig._generateId('ralinkpop');
            popnode.attr('id', popid);
            popper.ncbipopper({
                destSelector: "#" + popid,
                destPosition: 'top right', 
                triggerPosition: 'middle left', 
                hasArrow: true, 
                arrowDirection: 'right',
                isTriggerElementCloseClick: false,
                adjustFit: 'none',
                openAnimation: 'none',
                closeAnimation: 'none',
                delayTimeout : 130
            });
        }); // end each loop  
    });// end on page ready
})( jQuery );

Portal.Portlet.HistoryDisplay = Portal.Portlet.NCBIPageSection.extend({

	init: function(path, name, notifier) {
		console.info("Created History Ad...");
		this.base(path, name, notifier);    
	},
	
	send: {
      'Cmd': null      
    },   
    
    receive: function(responseObject, userArgs) {  
         var cmd = userArgs.cmd;
         var rootNode = document.getElementById('HTDisplay'); 
         var ul = document.getElementById('activity');
         var resp = responseObject.responseText;
             
         if (cmd == 'HTOn') { 
            rootNode.className = '';    // hide all msg and the turnOn link
            try {
            //alert(resp);
                // Handle timeouts
                if (responseObject.status == 408) { 
                    rootNode.className = 'HTOn'; // so that the following msg will show up
                    rootNode.innerHTML = "<p class='HTOn'>Your browsing activity is temporarily unavailable.</p>";
                    return;
                }
                   
                 // Looks like we got something...
                 resp = '(' + resp + ')';
                 var JSONobj = eval(resp);
                 
                 // Build new content (ul)
                 var newHTML = JSONobj.Activity;
                 var newContent = document.createElement('div');
                 newContent.innerHTML = newHTML;
                 var newUL = newContent.getElementsByTagName('ul')[0];
                 //alert(newHTML);
                 //alert(newContent.innerHTML);
                 //alert(newUL.innerHTML);
                 // Update content
                 rootNode.replaceChild(newUL, ul);
                 //XHR returns no activity (empty ul), e.g. activity cleared
                 if (newUL.className == 'hide')                     
                     rootNode.className = 'HTOn';  // show "Your browsing activity is empty." message
                 
            }         
            catch (e) {
                //alert('error');
                rootNode.className = 'HTOn'; // so that the following msg will show up
                rootNode.innerHTML = "<p class='HTOn'>Your browsing activity is temporarily unavailable.</p>";
           }
         }
         else if (cmd == 'HTOff') {                         
             if (ul != null) { 
                 ul.className='hide'; 
                 ul.innerHTML = ''; // clear activity
             }
             rootNode.className = 'HTOff';    // make "Activity recording is turned off." and the turnOn link show up             
         }
         else if (cmd == 'ClearHT') { 
             var goAhead = confirm('Are you sure you want to delete all your saved Recent Activity?');
             if (goAhead == true) { 
                 if ( rootNode.className == '') { //                 
                     rootNode.className = 'HTOn';  // show "Your browsing activity is empty." message                                  
                     if (ul != null) {
                         ul.className='hide'; 
                         ul.innerHTML = '';
                     }
                 }
             }
         } 
         
    },
    
	listen: {
	  'Cmd' : function(sMessage, oData, sSrc){
			console.info("Inside Cmd in HistoryDisplay: " + oData.cmd);
			this.setValue("Cmd", oData.cmd);
	  },	  
		
      "HistoryToggle<click>" : function(e, target, name){
         //alert(target.getAttribute("cmd"));
         this.send.Cmd({'cmd': target.getAttribute("cmd")});         
         console.info("Inside HistoryToggle in HistoryDisplay: " + target.getAttribute("cmd"));
         
         //var site = document.forms[0]['p$st'].value;
         var cmd =  target.getAttribute("cmd");     
               
         // Issue asynchronous call to XHR service, callback is to update the portlet output            
         this.doRemoteAction(target.getAttribute("cmd"));                      
      }, 
      
      "HistoryOn<click>" : function(e, target, name){
         this.send.Cmd({'cmd': target.getAttribute("cmd")});
         //$PN('Pubmed_ResultsSearchController').getInput('RecordingHistory').value = 'yes';		 
         console.info("Inside HistoryOn in HistoryDisplay: " + target.getAttribute("cmd"));
         this.doRemoteAction(target.getAttribute("cmd"));         
      },
      
      "ClearHistory<click>" : function(e, target, name){
         this.send.Cmd({'cmd': target.getAttribute("cmd")});
         this.doRemoteAction(target.getAttribute("cmd"));         
      }
    },
    
    'getPortletPath': function(){
        return this.realname + ".NCBIPageSection";
    }, 
    
    'doRemoteAction': function(command) {
         var site = document.forms[0]['p$st'].value;          
	     var resp = xmlHttpCall(site, this.realname, command, {}, this.receive, {'cmd': command}, this);
    }
});

;
Portal.Portlet.DbConnector = Portal.Portlet.extend({

	init: function(path, name, notifier) {
		var oThis = this;
		console.info("Created DbConnector");
		this.base(path, name, notifier);
		
		// reset Db value to original value on page load. Since LastDb is the same value as Db on page load and LastDb is not changed on
		// the client, this value can be used to reset Db. This is a fix for back button use.
		if (this.getValue("Db") != this.getValue("LastDb")){
		    this.setValue("Db", this.getValue("LastDb"));
		}
     
		// the SelectedIdList and id count from previous iteration (use a different attribute from IdsFromResult to prevent back button issues)
		Portal.Portlet.DbConnector.originalIdList = this.getValue("LastIdsFromResult");
		console.info("originalIdList " + Portal.Portlet.DbConnector.originalIdList);
		// if there is an IdList from last iteration set the count
		if (Portal.Portlet.DbConnector.originalIdList != ''){
			Portal.Portlet.DbConnector.originalCount = Portal.Portlet.DbConnector.originalIdList.split(/,/).length;
		}

		notifier.setListener(this, 'HistoryCmd', 
        	function(oListener, custom_data, sMessage, oNotifierObj) {
           		var sbTabCmd = $N(oThis.path + '.TabCmd');
           		sbTabCmd[0].value = custom_data.tab;
        	}
    		, null);
    
	},

	send: {
   		'SelectedItemCountChanged': null,
   		'newUidSelectionList': null,
   		'SavedSelectedItemCount': null,
   		'SavedUidList': null
	},

	listen: {
	
		//message from Display bar on Presentation change 
		'PresentationChange' : function(sMessage, oData, sSrc){
			
			// set link information only if it exists
			if (oData.dbfrom){
				console.info("Inside PresentationChange in DbConnector: " + oData.readablename);
				this.setValue("Db", oData.dbto);
				this.setValue("LinkSrcDb", oData.dbfrom);
				this.setValue("LinkName", oData.linkname);
				this.setValue("LinkReadableName", oData.readablename);
			}
			//document.forms[0].submit();
		},
		
		// various commands associated with clicking different form control elements
		'Cmd' : function(sMessage, oData, sSrc){
			console.info("Inside Cmd in DbConnector: " + oData.cmd);
			this.setValue("Cmd", oData.cmd);
			
			// back button fix, clear TabCmd
			if (oData.cmd == 'Go' || oData.cmd == 'PageChanged' || oData.cmd == 'FilterChanged' || 
			oData.cmd == 'DisplayChanged' || oData.cmd == 'HistorySearch' || oData.cmd == 'Text' || 
			oData.cmd == 'File' || oData.cmd == 'Printer' || oData.cmd == 'Order' || 
			oData.cmd == 'Add to Clipboard' || oData.cmd == 'Remove from Clipboard' || 
			oData.cmd.toLowerCase().match('details')){
				this.setValue("TabCmd", '');
				console.info("Inside Cmd in DbConnector, reset TabCmd: " + this.getValue('TabCmd'));
			}

		},
		
		
		// the term to be shown in the search bar, and used from searching
		'Term' : function(sMessage, oData, sSrc){
			console.info("Inside Term in DbConnector: " + oData.term);
			this.setValue("Term", oData.term);
		},
		
		
		// to indicate the Command Tab to be in
		'TabCmd' : function(sMessage, oData, sSrc){
			console.info("Inside TABCMD in DbConnector: " + oData.tab);
			this.setValue("TabCmd", oData.tab);
			console.info("DbConnector TabCmd: " + this.getValue("TabCmd"));
		},
		
		
		// message sent from SearchBar when db is changed while in a Command Tab
		'DbChanged' : function(sMessage, oData, sSrc){
			console.info("Inside DbChanged in DbConnector");
			this.setValue("Db", oData.db);
		},
		
		// Handles item select/deselect events
		// Argument is { 'id': item-id, 'selected': true or false }
		'ItemSelectionChanged' : function(sMessage, oData, oSrc) {
			var sSelection = this.getValue("IdsFromResult");
			var bAlreadySelected = (new RegExp("\\b" + oData.id + "\\b").exec(sSelection) != null);
	       	var count =0;
	       	
			if (oData.selected && !bAlreadySelected) {
				sSelection += ((sSelection > "") ? "," : "") + oData.id;
			   	this.setValue("IdsFromResult", sSelection);
			   	if (sSelection.length > 0){
			   		count = sSelection.split(',').length;
			   	}
			   	this.send.SelectedItemCountChanged({'count': count});
			   	this.send.newUidSelectionList({'list': sSelection});
			   	jQuery(document).trigger("itemsel",{'list': sSelection});
		   	} else if (!oData.selected && bAlreadySelected) {
				sSelection = sSelection.replace(new RegExp("^"+oData.id+"\\b,?|,?\\b"+oData.id+"\\b"), '');
		   	   	this.setValue("IdsFromResult", sSelection);
				console.info("Message ItemSelectionChanged - IdsFromResult after change:  " + this.getValue("IdsFromResult"));
			   	if (sSelection.length > 0){
			   		count = sSelection.split(',').length;
			   	}
				console.info("Message ItemSelectionChanged - IdsFromResult length:  " + count);   
				this.send.SelectedItemCountChanged({'count': count});
			   	this.send.newUidSelectionList({'list': sSelection});
			   	jQuery(document).trigger("itemsel",{'list': sSelection});
		   	}
		},
				
		// FIXME: This is the "old message" that is being phased out.
		// when result citations are selected, the list of selected ids are intercepted here,
		// and notification sent that selected item count has changed.
		'newSelection' : function(sMessage, oData, sSrc){
		
			// Check if we already have such IDs in the list
			var newList = new Array();
			var haveNow = new Array();
			if(Portal.Portlet.DbConnector.originalIdList){
				haveNow = Portal.Portlet.DbConnector.originalIdList.split(',');
				newList = haveNow;
			}
			
			var cameNew = new Array();
			if (oData.selectionList.length > 0) {
				cameNew = oData.selectionList;
			}
			
			if (cameNew.length > 0) {
				for(var ind=0;ind<cameNew.length;ind++) {
					var found = 0;
					for(var i=0;i<haveNow.length;i++) {
						if (cameNew[ind] == haveNow[i]) {
							found = 1;
							break;
						}
					}
						//Add this ID if it is not in the list
					if (found == 0) {
						newList.push(cameNew[ind]);
					}
				}
			}
			else {
				newList = haveNow;
			}

				// if there was an IdList from last iteration add new values to old
			var count = 0;
			if ((newList.length > 0) && (newList[0].length > 0)){
				count = newList.length;
			}
			
			console.info("id count = " + count);
			this.setValue("IdsFromResult", newList.join(","));
			
			this.send.SelectedItemCountChanged({'count': count});
			this.send.newUidSelectionList({'list': newList.join(",")});
			jQuery(document).trigger("itemsel",{'list': newList.join(",")});
		},


		// empty local idlist when list was being collected for other purposes.
		//used by Mesh and Journals (empty UidList should not be distributed, otherwise Journals breaks)
		// now used by all reports for remove from clipboard function.
		'ClearIdList' : function(sMessage, oData, sSrc){
			this.setValue("IdsFromResult", '');
			this.send.SelectedItemCountChanged({'count': '0'});
			this.send.newUidSelectionList({'list': ''});
			jQuery(document).trigger("itemsel",{'list': ""});
		}, 


		// back button fix: when search backend click go or hot enter on term field,
		//it also sends db. this db should be same as dbconnector's db
		'SearchBarSearch' : function(sMessage, oData, sSrc){
			if (this.getValue("Db") != oData.db){
				this.setValue("Db", oData.db);
			}
		},
		
		// back button fix: whrn links is selected from DisplayBar,
		//ResultsSearchController sends the LastQueryKey from the results on the page
		// (should not be needed by Entrez 3 code)
		'LastQueryKey' : function(sMessage, oData, sSrc){
			if (this.getInput("LastQueryKey")){
				this.setValue("LastQueryKey", oData.qk);
			}
		},
		
		'QueryKey' : function(sMessage, oData, sSrc){
			if (this.getInput("QueryKey")){
				this.setValue("QueryKey", oData.qk);
			}
		},
		
		
		//ResultsSearchController asks for the initial item count in case of send to file 
		'needSavedSelectedItemCount' : function(sMessage, oData, sSrc){
			var count = 0;
			if(this.getInput("IdsFromResult")){
				if (this.getValue("IdsFromResult").length > 0){
					count = this.getValue("IdsFromResult").split(',').length;
				}
				console.info("sending SavedSelectedItemCount from IdsFromResult: " + count);
			}
			else{
				count = Portal.Portlet.DbConnector.originalCount;
				console.info("sending SavedSelectedItemCount from OriginalCount: " + count);
			}
			this.send.SavedSelectedItemCount({'count': count});
		},
		
		// Force form submit, optionally passing db, term and cmd parameters
		'ForceSubmit': function (sMessage, oData, sSrc)
		{
		    if (oData.db)
    			this.setValue("Db", oData.db);
		    if (oData.cmd)
    			this.setValue("Cmd", oData.cmd);
		    if (oData.term)
    			this.setValue("Term", oData.term);
    		Portal.requestSubmit ();
		},
		
		'LinkName': function (sMessage, oData, sSrc){
		    this.setValue("LinkName", oData.linkname);
		},
		
		'IdsFromResult': function (sMessage, oData, sSrc){
		    this.setValue("IdsFromResult", oData.IdsFromResult);
		},
		
		'SendSavedUidList': function (sMessage, oData, sSrc){
		    this.send.SavedUidList({'idlist': this.getValue("IdsFromResult")});
		}
		
	}, //listen
	
	/* other portlet functions */
	
	// DisplayBar in new design wants selected item count
	'SelectedItemCount': function(){
	    var count = 0;
		if(this.getInput("IdsFromResult")){
			if (this.getValue("IdsFromResult") != ''){
				count = this.getValue("IdsFromResult").split(',').length;
			}
		}
		else{
			count = Portal.Portlet.DbConnector.originalCount;
		}
		return count;
	},
	
	'SelectedItemList': function(){
		if(this.getInput("IdsFromResult") && this.getValue("IdsFromResult") != ''){
			return this.getValue("IdsFromResult");
		}
		else{
			return Portal.Portlet.DbConnector.originalIdList;
		}
		
	},
	setValue: function(name, value){
	    if(name == 'Term')
	        value = jQuery.trim(value);
	    this.base(name,value);
	}
},
{
	originalIdList: '',
	originalCount: 0
});

function getEntrezSelectedItemCount() {
    return $PN('DbConnector').SelectedItemCount();
}

function getEntrezSelectedItemList() {
    return $PN('DbConnector').SelectedItemList();
}
