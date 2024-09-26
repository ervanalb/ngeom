console.debug("Loading wasm…");

import("./pkg").then(on_wasm_loaded).catch(on_error);

function on_wasm_loaded(pkg) {
    console.debug("Wasm loaded. Starting app…");

    let handle = new pkg.WebHandle();

    function check_for_panic() {
        if (handle.has_panicked()) {
            console.error("The app has crashed");

            // The demo app already logs the panic message and callstack, but you
            // can access them like this if you want to show them in the html:
            // console.error(`${handle.panic_message()}`);
            // console.error(`${handle.panic_callstack()}`);

            document.getElementById("physics_demo_canvas").remove();
            document.getElementById("center_text").innerHTML = `
                <p>
                    The app has crashed.
                </p>
                <p style="font-size:10px" align="left">
                    ${handle.panic_message()}
                </p>
                <p style="font-size:14px">
                    See the console for details.
                </p>
                <p style="font-size:14px">
                    Reload the page to try again.
                </p>`;
        } else {
            let delay_ms = 1000;
            setTimeout(check_for_panic, delay_ms);
        }
    }

    check_for_panic();

    handle.start(document.getElementById("physics_demo_canvas"))
      .then(on_app_started).catch(on_error);
}

function on_app_started(handle) {
    // Call `handle.destroy()` to stop. Uncomment to quick result:
    // setTimeout(() => { handle.destroy(); handle.free()) }, 2000)

    console.debug("App started.");
    document.getElementById("center_text").innerHTML = '';

    // Make sure the canvas is focused so it can receive keyboard events right away:
    document.getElementById("the_canvas_id").focus();
}

function on_error(error) {
    console.error("Failed to start: " + error);
    document.getElementById("the_canvas_id").remove();
    document.getElementById("center_text").innerHTML = `
        <p>
            An error occurred during loading:
        </p>
        <p style="font-family:Courier New">
            ${error}
        </p>
        <p style="font-size:14px">
            Make sure you use a modern browser with WebGL and WASM enabled.
        </p>`;
}
