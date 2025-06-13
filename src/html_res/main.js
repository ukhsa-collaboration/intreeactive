'use strict'

function showLabels() {
    isLabelsShown = true

    var labelLeftPadding = "\t\t\t"

    var labelArray

    if (builtLabelArray.length) {
        labelArray = targetElm
            .data[0]
            .text
            .map(x => x ? inputMetadata[x.split("<br")[0]] : x)
            .map(x => x ? `${labelLeftPadding}${builtLabelArray.map(labelField => x[labelField]).join(builtLabelSeparator)}` : "")
    }
    else {
        labelArray = targetElm
            .data[0]
            .text
            .map(x => x ? inputMetadata[x.split("<br")[0]] : x)
            .map(x => x ? `${labelLeftPadding}${x[labelBy]}` : "")
    }

    // we need to nest labelArray in another Array due to
    // a Plotly WebGL bug where it only reads the first
    // element of the input array.
    // The non-WebGL version does not need this nesting.
    // We could check if the plot is of type ScatterGl...
    var updateData = {
        text: [labelArray],
        visible: true
    }

    if (rainbowMode) {
        updateData["textfont"] = {
            color: targetElm.data[0].marker.color
        }
    }

    Plotly.restyle(targetElm, updateData, 1);
}

function hideLabels() {
    isLabelsShown = false

    Plotly.restyle(targetElm, { "visible": false }, 1);
}

function doHighlight(inputString) {
    var colArray = targetElm.data[0].marker.color
    var labelArray = targetElm.data[0].text.map(x => x ? x.split("<br>")[0] : x)
    var data = {}

    // store the original colours
    // if we haven't already changed them
    if (!isCustomColoursEnabled) {
        isCustomColoursEnabled = true
        originalColours = colArray
    }

    if (inputString.length) {
        data = {
            marker: {
                color: labelArray.map(x =>
                    x
                        ? x.includes(inputString)
                            ? "red"
                            : "rgb(100, 100, 100)"
                        : "rgb(100,100,100)"
                ),
                size: labelArray.map(x => x ? x.includes(inputString) ? 15 : 10 : x)
            }
        }
    }
    else {
        // data = {
        //     marker: {
        //         color: labelArray.map(x => "rgb(100, 100, 100)")
        //     }
        // }
        restoreOriginalColours()
        return
    }

    Plotly.restyle(targetElm, data);
}

function doHighlightIdArray(inputArray, inputId) {
    // exact matches to valid IDs only
    // inputArray gets colour one
    // inputId if present gets colour two
    var colArray = targetElm.data[0].marker.color
    var labelArray = targetElm.data[0].text.map(x => x ? x.split("<br>")[0] : x)
    var data = {}
    var colourOne = "blue"
    var colourTwo = "red"

    // store the original colours
    // if we haven't already changed them
    if (!isCustomColoursEnabled) {
        isCustomColoursEnabled = true
        originalColours = colArray
    }

    if (inputArray.length) {
        data = {
            marker: {
                color: labelArray.map(x => x ? x == inputId ? colourTwo : inputArray.includes(x) ? colourOne : "rgb(100, 100, 100)" : "rgb(100, 100, 100)"),
                size: labelArray.map(x => x ? x.includes(inputId) ? 15 : 10 : x)
            }
        }
    }
    else {
        // data = {
        //     marker: {
        //         color: labelArray.map(x => "rgb(100, 100, 100)")
        //     }
        // }
        restoreOriginalColours()
        return
    }

    Plotly.restyle(targetElm, data);
}

function getCurrentDropdownStatus() {
    var activeColourIdx = targetElm.layout.updatemenus[0].active
    var activeColourCategory = targetElm.layout.updatemenus[0].buttons[activeColourIdx]

    return activeColourCategory
}


function restoreOriginalColours() {
    if (!originalColours) {
        console.warn("nothing to restore")
        return
    }

    var data = {
        marker: {
            color: originalColours,
            size: originalMarkerSize
        }
    }

    isCustomColoursEnabled = false

    Plotly.restyle(targetElm, data);
}

function populateMetadataTableByIdPair(inputId, inputId2) {
    var tempMetadata = inputMetadata[inputId]
    var tempMetadata2 = inputMetadata[inputId2]

    // pre-emptively blank the div
    metadataElm.innerHTML = ""


    if (tempMetadata && tempMetadata2) {
        // if we have both, do a paired table
        var tempMetadataElm = formatPairedKeyValuesAsNiceTable(tempMetadata, tempMetadata2)
        tempMetadataElm.classList.add("metadataTable")
        metadataElm.append(tempMetadataElm)
    }
    else if (tempMetadata || tempMetadata2) {
        // if we only have one, do a table with just that one
        var tempMetadataElm = formatKeyValuesAsNiceTable(tempMetadata ? tempMetadata : tempMetadata2)
        tempMetadataElm.classList.add("metadataTable")
        metadataElm.append(tempMetadataElm)
    }
    else {
        // there's nothing to do
        console.error("Failed to grab metadata for inputId", inputId)
    }
}

class NearestNeighbourID {
    constructor(inputID) {
        this.nnID = inputID
        this.elm = document.createElement("div")

        this.elm.innerHTML = inputID
        this.elm.classList.add("nnID")
        this.elm.dataset.nodeId = inputID
        this.elm.setAttribute("onclick", "doHighlight(this.dataset.nodeId)")
    }
}

class NearestNeighbourSNPCount {
    constructor(inputCount) {
        this.nnSNPs = inputCount
        this.elm = document.createElement("div")
        this.isPlural = this.nnSNPs != 1

        this.elm.innerHTML = `${this.nnSNPs} SNP${this.isPlural ? "s" : ""}`
        this.elm.classList.add("nnSNPs")
    }
}

class NearestNeighbourCount {
    constructor(inputCount) {
        this.nnCount = inputCount
        this.elm = document.createElement("div")
        this.isPlural = this.nnCount != 1

        this.elm.innerHTML = `${this.nnCount} Neighbour${this.isPlural ? "s" : ""}`
        this.elm.classList.add("nnCount")
    }
}

class NearestNeighbourWithMeta {
    constructor(inputId, inputSnps) {
        this.nnId = inputId
        this.nnCount = inputSnps
        this.elm = document.createElement("div")
        this.isPlural = this.nnSNPs != 1

        var nnIDelm = document.createElement("div")
        nnIDelm.classList.add("nnID")
        nnIDelm.innerHTML = this.nnId
        nnIDelm.dataset.nodeId = inputId
        nnIDelm.setAttribute("onclick", "doHighlight(this.dataset.nodeId)")

        var nnSNPselm = document.createElement("div")
        nnSNPselm.classList.add("nnSNPs")
        nnSNPselm.innerHTML = this.nnCount

        this.elm.classList.add("nnWithMeta")
        this.elm.append(nnIDelm, nnSNPselm)
    }
}

class HighlightAllBtn {
    constructor() {
        this.elm = document.createElement("div")

        this.elm.classList.add("nnHighlightAll")
        this.elm.innerHTML = "Highlight all"
        this.elm.setAttribute("onclick", "doHighlightIdArray(this.idArray, this.mainId)")
        this.elm.idArray = []
    }
}

function formatKeyValuesAsNiceTable(inputObj) {
    var tempTable = document.createElement("table")
    var tempTableBody = document.createElement("tbody")

    for (var k in inputObj) {
        var tempRow = document.createElement("tr")
        var th_key = document.createElement("th")
        var td_value = document.createElement("td")

        th_key.innerHTML = k.replaceAll("_", " ")
        td_value.innerHTML = inputObj[k]

        tempRow.append(th_key, td_value)
        tempTableBody.append(tempRow)

        if (k == "ID") {
            // add a class that we can use for styling
            td_value.classList.add("metadataId")

            // store the ID in the node
            td_value.dataset.nodeId = inputObj[k]

            // add the highlight onclick event
            td_value.setAttribute("onclick", `doHighlight(this.dataset.nodeId)`)
        }

        if (k == "Nearest_neighbour") {
            // overwrite values
            td_value.innerHTML = ""

            // set inline CSS on row
            td_value.parentElement.style.verticalAlign = "top"

            // stage highlight all button
            var haBtn = new HighlightAllBtn()
            haBtn.elm.mainId = inputObj.ID

            if (!snpThreshold) {
                // it's an array, so we want to treat it differently
                var tempVal = inputObj[k]
                var snpCount = parseInt(tempVal[0].split("=").slice(-1))

                td_value.append(
                    new NearestNeighbourCount(tempVal.length).elm,
                    // new NearestNeighbourSNPCount(snpCount).elm
                )

                var tempRecords = tempVal.map(x => x.split("="))

                tempRecords.map(x => new NearestNeighbourWithMeta(x[0], x[1])).forEach(x => {
                    td_value.append(x.elm)
                })

                tempRecords.forEach(x => { haBtn.elm.idArray.push(x[0]) })

                if (tempRecords.length > 1) {
                    td_value.append(haBtn.elm)
                }
            }
            // if we've got a defined SNP threshold
            else {
                var tempId = inputObj["ID"]
                var nnTempId = getNeighboursWithinSnpThreshold(tempId, snpThreshold)

                if (nnTempId.length) {
                    td_value.append(...nnTempId.map(x => new NearestNeighbourWithMeta(x[0], x[1]).elm))
                    haBtn.elm.idArray.push(...nnTempId.map(x => x[0]))
                    if (nnTempId.length > 1) {
                        td_value.append(haBtn.elm)
                    }
                }
                else {
                    td_value.innerHTML = "No neighbours within threshold"
                }
            }
        }
    }

    tempTable.append(tempTableBody)

    return tempTable
}

function formatPairedKeyValuesAsNiceTable(inputObj, inputObj2) {
    var tempTable = document.createElement("table")
    var tempTableHead = document.createElement("thead")
    var tempTableBody = document.createElement("tbody")

    for (var k in inputObj) {
        var tempRow = document.createElement("tr")
        var th_key = document.createElement("th")
        var td_value = document.createElement("td")
        var td_value2 = document.createElement("td")

        th_key.innerHTML = k.replaceAll("_", " ")
        td_value.innerHTML = inputObj[k]
        td_value2.innerHTML = inputObj2[k]

        tempRow.append(th_key, td_value, td_value2)
        tempTableBody.append(tempRow)

        if (k == "ID") {
            // add a class that we can use for styling
            td_value.classList.add("metadataId")
            td_value2.classList.add("metadataId")

            // store the ID in the node
            td_value.dataset.nodeId = inputObj[k]
            td_value2.dataset.nodeId = inputObj2[k]

            // add the highlight onclick event
            td_value.setAttribute("onclick", `doHighlight(this.dataset.nodeId)`)
            td_value2.setAttribute("onclick", `doHighlight(this.dataset.nodeId)`)

            // add the snp distances row
            var snpDistRow = `<tr><th>SNP distance</th><td colspan=2>${getSnpDistance(inputObj[k], inputObj2[k])}</td></tr>`
            tempTableBody.innerHTML += snpDistRow
        }

        if (k == "Nearest_neighbour") {
            // overwrite values
            td_value.innerHTML = ""
            td_value2.innerHTML = ""

            // set inline CSS on row
            td_value.parentElement.style.verticalAlign = "top"

            // stage highlight all button
            var haBtn = new HighlightAllBtn()
            haBtn.elm.mainId = inputObj.ID
            var haBtn2 = new HighlightAllBtn()
            haBtn2.elm.mainId = inputObj2.ID

            // if we're just using the minimum SNP distance
            // as our nearest neighbours
            if (!snpThreshold) {
                // it's an array, so we want to treat it differently
                var tempVal = inputObj[k]
                var tempVal2 = inputObj2[k]
                var snpCount = parseInt(tempVal[0].split("=").slice(-1))
                var snpCount2 = parseInt(tempVal2[0].split("=").slice(-1))

                td_value.append(
                    new NearestNeighbourCount(tempVal.length).elm,
                    // new NearestNeighbourSNPCount(snpCount).elm
                )

                var tempRecords = tempVal.map(x => x.split("="))

                tempRecords.map(x => new NearestNeighbourWithMeta(x[0], x[1])).forEach(x => {
                    td_value.append(x.elm)
                })

                tempRecords.forEach(x => { haBtn.elm.idArray.push(x[0]) })

                if (tempRecords.length > 1) {
                    td_value.append(haBtn.elm)
                }

                td_value2.append(
                    new NearestNeighbourCount(tempVal2.length).elm,
                    // new NearestNeighbourSNPCount(snpCount2).elm
                )

                var tempRecords2 = tempVal2.map(x => x.split("="))

                tempRecords2.map(x => new NearestNeighbourWithMeta(x[0], x[1])).forEach(x => {
                    td_value2.append(x.elm)
                })

                tempRecords2.forEach(x => { haBtn2.elm.idArray.push(x[0]) })
            }
            // if we've got a defined SNP threshold
            else {
                var tempId = inputObj["ID"]
                var tempId2 = inputObj2["ID"]
                var nnTempId = getNeighboursWithinSnpThreshold(tempId, snpThreshold)
                var nnTempId2 = getNeighboursWithinSnpThreshold(tempId2, snpThreshold)

                if (nnTempId.length) {
                    td_value.append(...nnTempId.map(x => new NearestNeighbourWithMeta(x[0], x[1]).elm))
                    haBtn.elm.idArray.push(...nnTempId.map(x => x[0]))
                    if (nnTempId.length > 1) {
                        td_value.append(haBtn.elm)
                    }
                }
                else {
                    td_value.innerHTML = "No neighbours within threshold"
                }

                if (nnTempId2.length) {
                    td_value2.append(...nnTempId2.map(x => new NearestNeighbourWithMeta(x[0], x[1]).elm))
                    haBtn2.elm.idArray.push(...nnTempId2.map(x => x[0]))
                    if (nnTempId2.length > 1) {
                        td_value2.append(haBtn2.elm)
                    }
                }
                else {
                    td_value2.innerHTML = "No neighbours within threshold"
                }
            }
        }
    }

    tempTable.append(tempTableBody)

    return tempTable
}

function getNeighboursWithinSnpThreshold(inputId, inputMaxSnps) {
    var idx = inputSnpMatrix.index.indexOf(inputId)
    var data = inputSnpMatrix.data[idx]
    var returnIds = []

    data.forEach((x, i) => {
        // if below or equal to the threshold, and
        // not the input ID
        if ((x <= inputMaxSnps) && (i != idx)) {
            // push the ID and the number of SNPs
            returnIds.push([inputSnpMatrix.index[i], x])
        }
    })

    // sort by SNP count
    // doesn't sort by label
    returnIds.sort(function (a, b) { return a[1] - b[1] })

    return returnIds
}

function getParentTableElm(inputElm) {
    try {
        var parent = inputElm.parentElement
        if (parent.tagName == "TABLE") {
            return parent
        }
        else {
            return getParentTableElm(parent)
        }
    }
    catch {
        console.error("Failed to find table parent for inputElm", inputElm)
        return
    }
}

// function doHighlightIdAndNeighboursinColumn(inputElm){
//     var parentTableElm = getParentTableElm(inputElm)
//     var colIdx = inputElm.parentElement
// }

function getSnpDistance(inputId1, inputId2) {
    if (!inputId1 || !inputId2) {
        console.warn("Two ids required to get snp distance. Input:", inputId1, inputId2)
        return
    }

    try {
        var idx1 = inputSnpMatrix.index.indexOf(inputId1)
        var idx2 = inputSnpMatrix.index.indexOf(inputId2)

        return inputSnpMatrix.data[idx1][idx2]
    }
    catch {
        console.warn("Two valid ids required to get snp distance. Input:", inputId1, inputId2)
        return
    }
}

function initIdDatalist() {
    var inputElm1 = document.createElement("input")
    var inputElm2 = document.createElement("input")
    var dataListElm = document.createElement("datalist")
    var label1 = document.createElement("label")
    var label2 = document.createElement("label")
    var clearBtnElm = document.createElement("button")

    dataListElm.id = "metaIdsDatalist";
    inputElm1.setAttribute("list", "metaIdsDatalist");
    inputElm1.id = "metaIdsDatalistInput1";
    inputElm1.setAttribute("oninput", "populateMetadataTableByIdPair(this.value, document.getElementById('metaIdsDatalistInput2').value)");
    inputElm1.setAttribute("type", "text")
    inputElm2.setAttribute("list", "metaIdsDatalist");
    inputElm2.setAttribute("type", "text")
    inputElm2.id = "metaIdsDatalistInput2";
    inputElm2.setAttribute("oninput", "populateMetadataTableByIdPair(document.getElementById('metaIdsDatalistInput1').value, this.value)");
    label1.setAttribute("for", "metaIdsDatalistInput1");
    label2.setAttribute("for", "metaIdsDatalistInput2");

    // label1.innerHTML = "IDs for comparison";
    // label2.innerHTML = "Select ID 2:";

    inputElm1.setAttribute("placeholder", "ID 1");
    inputElm2.setAttribute("placeholder", "ID 2");

    clearBtnElm.innerHTML = "Clear";
    clearBtnElm.setAttribute("onclick", "clearMetadataIdsAndResults()");

    // don't forget to put a ; after the previous line
    // the spread syntax operator complains otherwise
    var tempArray = [...Object.keys(inputMetadata)]
    tempArray.sort() // toSorted not supported in old firefox
    tempArray.forEach(x => {
        var tempElm = document.createElement("option")
        tempElm.value = x
        dataListElm.appendChild(tempElm)
    })

    metadataControlsElm.append(label1, inputElm1, label2, inputElm2, dataListElm, clearBtnElm)
}

function clearMetadataIdsAndResults() {
    document.getElementById("metaIdsDatalistInput1").value = ""
    document.getElementById("metaIdsDatalistInput2").value = ""
    metadataElm.innerHTML = ""
    restoreOriginalColours()
}

function initLabelByDropdown() {
    var metadataFields = [...Object.keys(inputMetadata[Object.keys(inputMetadata)[0]])]
    var dropdownElm = document.createElement("select")
    var excludeFields = ["Nearest_neighbour"]

    dropdownElm.id = "labelByDropdown"
    dropdownElm.setAttribute("onchange", "labelByDropdownChange(this.value)")

    metadataFields.forEach(x => {
        if (!excludeFields.includes(x)) {
            var tempElm = document.createElement("option")
            tempElm.innerHTML = x
            tempElm.value = x
            dropdownElm.append(tempElm)
        }
    })

    var addFieldBtn = document.createElement("button")
    addFieldBtn.id = "addFieldBtn"
    addFieldBtn.innerHTML = "+"
    addFieldBtn.setAttribute("onclick", "addLabelField()")

    var builtLabelFieldDiv = document.createElement("div")
    builtLabelFieldDiv.id = "builtLabelField"

    document.getElementById("labelDropdownCell").append(dropdownElm, addFieldBtn, builtLabelFieldDiv)
}

function labelByDropdownChange(inputField) {
    // different behaviour if we're using the builder
    if (builtLabelArray.length) {
        return
    }

    labelBy = inputField
    if (isLabelsShown) {
        showLabels()
    }
}

function addLabelField() {
    // get current label dropdown value, add it to global array
    builtLabelArray.push(document.getElementById("labelByDropdown").value)

    // do update
    updateLabelField()
}

function updateLabelField() {
    var builtLabelFieldElm = document.getElementById("builtLabelField")

    // update builtLabelField div
    builtLabelFieldElm.innerHTML = "" // blank it
    builtLabelArray.forEach((x, i) => {
        var tempElm = document.createElement("div")
        tempElm.classList.add("labelBuilderLabel")
        tempElm.innerHTML = x
        tempElm.labelIdx = i
        tempElm.setAttribute("onclick", "removeLabelField(this)")
        builtLabelFieldElm.append(tempElm)
    })

    if (isLabelsShown) {
        showLabels()
    }
}

function removeLabelField(inputElm) {
    // get idx
    var idxOfElm = inputElm.labelIdx

    // remove element from array
    builtLabelArray.splice(idxOfElm)

    // remove element from DOM
    inputElm.remove()

    // don't use this method, it has a race condition of some kind
    // just use the remove from DOM method
    // // update
    // updateLabelField()

    if (isLabelsShown) {
        showLabels()
    }
}

function refreshMetadataTable() {
    var tempVal1 = document.getElementById('metaIdsDatalistInput1').value
    var tempVal2 = document.getElementById('metaIdsDatalistInput2').value

    if (tempVal1.length || tempVal2.length) {
        populateMetadataTableByIdPair(
            tempVal1,
            tempVal2
        )
    }
    else {
        console.log("Nothing to refresh")
    }
}

function snpThresholdChange(inputValue) {
    switch (inputValue) {
        case "Use min":
            snpThreshold = false
            refreshMetadataTable()
            break
        case "Set manually":
            snpThreshold = parseInt(snpThresholdSpinnerElm.value)
            refreshMetadataTable()
            break
        default:
            console.error("Error in snpThresholdChange, unrecognised option", inputValue)
    }
}

function snpThresholdSpinnerChange(inputValue) {
    document.querySelector("input[value='Set manually']").checked = true
    snpThresholdChange('Set manually')
}

function initSnpThresholdRadio() {
    var radioElms = [...document.querySelectorAll('input[name="snpThresholdRadio"]')]

    radioElms.forEach(x => {
        x.setAttribute("onchange", "snpThresholdChange(this.value)")
    })

}

function initHelpToggle() {
    var toggleElm = document.createElement("div")
    toggleElm.id = "helpToggle"

    toggleElm.innerHTML = "?"
    toggleElm.setAttribute("onclick", "toggleHelp()")

    document.body.append(toggleElm)
}

function toggleHelp() {
    document.getElementById("helpModal").classList.toggle("hidden")
}

function initDarkModeToggle() {
    var toggleElm = document.createElement("div")
    toggleElm.id = "darkModeToggle"

    toggleElm.innerHTML = sunEmoji
    toggleElm.setAttribute("onclick", "toggleDarkMode()")

    document.body.append(toggleElm)
}

function toggleDarkMode() {
    isDarkModeEnabled = !isDarkModeEnabled

    if (isDarkModeEnabled && !document.querySelectorAll("#darkModeStyle").length) {
        var tempElm = document.createElement("style")
        tempElm.id = "darkModeStyle"
        tempElm.innerHTML = `
            :root {
                --page-background: #353144;
                --panel-background : #433e56;
                --panel-foreground: white;
                --table-row-odd: #433e56;
                --table-row-even: #716799;
                --table-row-text: white;
                --main-border: #595565;
            }
            .plot-container {
                filter: invert(85%) hue-rotate(180deg);
            }
            #metadataDivControls{
                background: var(--table-row-even);
            }
            #darkModeToggle {
                background: var(--page-background);
                color: white;
                border: 5px solid var(--table-row-even);
                box-shadow: 0px 0px 4px var(--panel-background);
            }
            #helpToggle {
                background: var(--page-background);
                color: white;
                border: 5px solid var(--table-row-even);
                box-shadow: 0px 0px 4px var(--panel-background);
            }
            #builtLabelField .labelBuilderLabel{
                background: white;
                color: black;
                border: 1px solid var(--table-row-odd);
            }
        `
        document.body.append(tempElm)
    }
    else {
        document.querySelectorAll("#darkModeStyle").forEach(x => x.remove())
    }
}

function setRangeX(inputMin, inputMax) {
    var layout = {
        xaxis: {
            range: [inputMin, inputMax],
            showgrid: false,
            zeroline: false
        }
    }

    Plotly.relayout(targetElm, layout);
}

function setRangeY(inputMin, inputMax) {
    var layout = {
        yaxis: {
            range: [inputMin, inputMax],
            showgrid: false,
            zeroline: false,
            visible: false
        }
    }

    Plotly.relayout(targetElm, layout);
}

function doZoomX(zoomOut) {
    var xRange = targetElm.layout.xaxis.range
    var xRangeDelta = xRange[1] - xRange[0]
    var scaleFactor = 0.1
    var newDelta = scaleFactor * xRangeDelta
    var newXmin = xRange[0]
    var newXmax = xRange[1]

    if (zoomOut) {
        newXmin -= newDelta
        newXmax += newDelta
    }
    else {
        newXmin += newDelta
        newXmax -= newDelta
    }

    setRangeX(newXmin, newXmax)
}

function doZoomY(zoomOut) {
    var yRange = targetElm.layout.yaxis.range
    var yRangeDelta = yRange[1] - yRange[0]
    var scaleFactor = 0.1
    var newDelta = scaleFactor * yRangeDelta
    var newYmin = yRange[0]
    var newYmax = yRange[1]

    if (zoomOut) {
        newYmin -= newDelta
        newYmax += newDelta
    }
    else {
        newYmin += newDelta
        newYmax -= newDelta
    }

    setRangeY(newYmin, newYmax)
}

function resetScale() {
    setRangeX(...initialAxisRange.xaxis)
    setRangeY(...initialAxisRange.yaxis)
}

function autoResizePlotHeight() {
    // detect available height and force
    // plot to relayout with that height
    var targetHeight = 0.9 * leftPanelElm.offsetHeight
    var minHeight = 300

    Plotly.relayout(targetElm, { height: Math.max(targetHeight, minHeight), width: leftPanelElm.offsetWidth })
}

function updatePlotTitle(newTitleString){
    var newTitleObj = {
        text: newTitleString,
        y: initialPlotTitleObj.y,
        yanchor: initialPlotTitleObj.yanchor
    }

    Plotly.relayout(targetElm, newTitleObj)
}

// globals
const targetElm = document.getElementsByClassName("plotly-graph-div")[0]
const metadataElm = document.getElementById("metadataDiv")
const metadataControlsElm = document.getElementById("metadataDivControls")
const highlightInputElm = document.getElementById("highlightInput")
const snpThresholdSpinnerElm = document.getElementById("snpThresholdSpinner")
const leftPanelElm = document.getElementById("leftPanel")
const initialAxisRange = {
    xaxis: targetElm.layout.xaxis.range,
    yaxis: targetElm.layout.yaxis.range
}
const initialPlotTitleObj = targetElm.layout.title
var originalColours
var originalMarkerSize = targetElm.data[0].marker.size
var isCustomColoursEnabled = false
var labelBy = "ID"
var isLabelsShown = false
var snpThreshold
const sunEmoji = "&#x2600;"
const moonEmoji = "&#x263E;"
var isDarkModeEnabled = false;
var rangeMin = -0.1
var rangeMax = 0.5
var maxInputElm = document.getElementById("maxInput")
var minInputElm = document.getElementById("minInput")
var rainbowMode = false
var builtLabelArray = []
var builtLabelSeparator = " | " // consider underscore, interpunct, space, slash...

// init
function init() {
    initIdDatalist();
    initLabelByDropdown();
    initSnpThresholdRadio();
    initDarkModeToggle();
    initHelpToggle();
    targetElm.on("plotly_update", function () {
        originalColours = targetElm.data[0].marker.color
    })
    autoResizePlotHeight();
    window.addEventListener("resize", autoResizePlotHeight);
}

init();