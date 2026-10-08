/* cellpacman Shiny application script.
 * Attached to the UI by app_ui() through an htmltools::htmlDependency().
 *
 * The server enables each analysis tab only after the stage it depends on
 * has completed (Data -> Dimension Reduction -> Clustering -> Trajectory).
 * It does so by sending a `setNavDisabled` custom message; this handler
 * toggles the Bootstrap `disabled` state on the matching navbar link. */
(function () {
  function findNavLink(navId, value) {
    return document.querySelector(
      '#' + navId + ' a.nav-link[data-value="' + value + '"]'
    );
  }

  function registerHandler() {
    if (!window.Shiny) {
      window.setTimeout(registerHandler, 50);
      return;
    }

    Shiny.addCustomMessageHandler('setNavDisabled', function (message) {
      var link = findNavLink(message.navId, message.value);
      if (!link) {
        return;
      }

      if (message.disabled) {
        link.classList.add('disabled');
        link.setAttribute('aria-disabled', 'true');
        link.setAttribute('tabindex', '-1');
      } else {
        link.classList.remove('disabled');
        link.removeAttribute('aria-disabled');
        link.removeAttribute('tabindex');
      }
    });
  }

  registerHandler();
})();
