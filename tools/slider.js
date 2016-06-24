el = document.getElementById("simple-ss");
el.onclick = links;

function links() {
  left = parseInt(getComputedStyle(el).getPropertyValue("background-position").split(" ", 1));

  /* DEFINE POSITIONS FOR CLICK EVENTS */
  if (left >= -450) {

    // First until about half way scrolled over
    //alert("first");
    //window.open("http://www.google.com");

  } else if (left >= -1350) {

    // Second when half way scrolled either side
    //alert("second");
    window.open("https://www.cnio.es/eventos/index.asp?ev=1&cev=136");

  } 

}