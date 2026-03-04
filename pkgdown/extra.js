/**
 * Tabsets for pkgdown articles (Bootstrap 5 compatible).
 * Adapted from rmarkdown tabsets.js for Pandoc .tabset and Bootstrap 5.
 */
(function() {
  "use strict";

  function buildTabsets(tocID) {
    // Find tabsets: Pandoc 2.8+ puts .tabset on heading; 2.7 and earlier on div
    var tabsets = document.querySelectorAll("div.section.tabset");
    if (tabsets.length === 0) {
      tabsets = Array.from(document.querySelectorAll("div.section")).filter(function(div) {
        return div.querySelector(":scope > .tabset");
      });
    }

    tabsets.forEach(function(tabsetEl) {
      var tabset = tabsetEl;
      var match = (tabset.className || "").match(/level(\d)/);
      if (!match) return;
      var tabsetLevel = parseInt(match[1], 10);
      var tabLevel = tabsetLevel + 1;

      var allTabs = tabset.querySelectorAll("div.section.level" + tabLevel);
      var tabs = Array.prototype.filter.call(allTabs, function(t) {
        var h = t.querySelector("h" + tabLevel);
        var isCode = t.querySelector(".sourceCode, pre.sourceCode");
        return h && !isCode;
      });
      if (!tabs.length) return;

      var navClass = tabset.classList.contains("tabset-pills") ? "nav-pills" : "nav-tabs";
      var tabList = document.createElement("ul");
      tabList.className = "nav " + navClass;
      tabList.setAttribute("role", "tablist");

      var tabContent = document.createElement("div");
      tabContent.className = "tab-content";

      var activeTab = 0;
      var firstTab = tabs[0];

      tabset.insertBefore(tabContent, firstTab);
      tabset.insertBefore(tabList, firstTab);

      for (var i = 0; i < tabs.length; i++) {
        var tab = tabs[i];
        var id = tab.id || ("tab-" + Math.random().toString(36).slice(2));
        if (tab.classList.contains("active")) activeTab = i;
        id = id.replace(/[.\/?&!#<>]/g, "").replace(/\s/g, "_");
        tab.id = id;

        var heading = tab.querySelector("h" + tabLevel);
        var headingText = heading ? heading.innerHTML : "Tab " + (i + 1);
        if (heading) heading.remove();

        var li = document.createElement("li");
        li.className = "nav-item";
        li.setAttribute("role", "presentation");
        var a = document.createElement("a");
        a.className = "nav-link";
        a.setAttribute("role", "tab");
        a.setAttribute("data-bs-toggle", "tab");
        a.setAttribute("href", "#" + id);
        a.setAttribute("aria-controls", id);
        a.innerHTML = headingText;
        li.appendChild(a);
        tabList.appendChild(li);

        tab.setAttribute("role", "tabpanel");
        tab.classList.add("tab-pane", "tabbed-pane");
        if (tabset.classList.contains("tabset-fade")) tab.classList.add("fade");
        tabContent.appendChild(tab);

        if (tocID && document.getElementById(tocID)) {
          var tocLink = document.querySelector("#" + tocID + " li a[href='#" + id + "']");
          if (tocLink && tocLink.parentNode) tocLink.parentNode.remove();
        }
      }

      var firstLi = tabList.children[activeTab];
      var firstPane = tabContent.children[activeTab];
      if (firstLi) firstLi.querySelector("a").classList.add("active");
      if (firstPane) {
        firstPane.classList.add("active", "show");
        if (tabset.classList.contains("tabset-fade")) firstPane.classList.add("in");
      }
    });
  }

  document.addEventListener("DOMContentLoaded", function() {
    buildTabsets("TOC");
  });
})();
