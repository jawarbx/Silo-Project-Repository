
.removeChild(document.querySelector('[id^="ack"]'))
list = arguments[0].querySelectorAll('[id^="p0"]')
text = []
for (let i=0; i<list.length; i++) {
    text.push(list[i].textContent)
}
return text.toString()