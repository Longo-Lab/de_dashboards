$(document).ready(function() {
  
  // Scrolling for sidebar
  $(window).scroll(function(){
    var h = 50,
        $scrollTop = $(window).scrollTop();

    $('.main-sidebar').css({
      'paddingTop': ($scrollTop > h ? 0 : h - $scrollTop) + 'px'
    });
  })
  
  // Buttons for modules enrichment tab
  $(document).on('click', '#biodomains a', function() {
      var $this = $(this),
          $value = $this.data('value');

      if ($this.hasClass('active')) return;
        
      $this.parent().find('a').removeClass('active');
      $this.addClass('active');
      
      $('#biodomains ~ img').attr('src', $value);
  })
  
  // Collapsible tabs
  $(document).on('click', '.nav-tabs .tab-group', function(e) {
      e.preventDefault();
      
      $(`.nav-tabs-custom > .nav-tabs .${$(this).text().trim().toLowerCase()}:not(.tab-group)`).animate({
            width: 'toggle',
            opacity: 'toggle'
        }, 500);
  })
  
});
