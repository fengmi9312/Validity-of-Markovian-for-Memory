
function openTab(tabName) {
    // Hide all tab content
    var tabContent = document.getElementsByClassName("tab-content");
    for (var i = 0; i < tabContent.length; i++) {
      tabContent[i].classList.remove("active");
    }
  
    // Show the selected tab content
    document.getElementById(tabName).classList.add("active");
  }
  
  