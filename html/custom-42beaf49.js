document.addEventListener("DOMContentLoaded", function () {
  var blockquotes = document.querySelectorAll("blockquote");
  blockquotes.forEach(function (bq) {
    var text = bq.innerText.trim().toUpperCase();
    if (text.startsWith("NOTE") || text.startsWith("!NOTE")) {
      bq.classList.add("alert-note");
    } else if (text.startsWith("WARNING") || text.startsWith("!WARNING")) {
      bq.classList.add("alert-warning");
    } else if (text.startsWith("IMPORTANT") || text.startsWith("!IMPORTANT")) {
      bq.classList.add("alert-important");
    } else if (text.startsWith("TIP") || text.startsWith("!TIP")) {
      bq.classList.add("alert-tip");
    }
  });
});
