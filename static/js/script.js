function changeFileName(){
    let name = document.getElementById('fileName');
    let f = document.cookie.split(';');
    for (let cookie of f) {
        let parts = cookie.trim().split('=');
        if (parts[0] === 'filename') {
            name.innerHTML = parts[1];
            break
        }
    }
}

function showGIF(){
    let non = document.getElementById('non-loading');
    non.style.opacity = "0.1";
    let loading = document.getElementById("loading");
    loading.style.display = "block";
}

function copy(id){
    let copyText = document.getElementById(id).textContent;
    if (copyText === ''){
        console.log('here!');
        copyText = document.getElementById(id).value;
    }
    let textArea = document.createElement('textarea');
    textArea.id = 'temp';
    textArea.textContent = copyText;
    document.body.append(textArea);
    textArea.select();
    document.execCommand("copy");
    /* Alert the copied text */
    alert("Copied to clipboard");
    document.getElementById('temp').remove()
}
