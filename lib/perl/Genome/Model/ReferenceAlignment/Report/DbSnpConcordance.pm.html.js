// hard code color indices to prevent from shifting as plots are toggled
var i = 0;
$.each(datasets, function(key, val) {
    val.color = i;
    ++i;
});

// insert checkboxes
var choiceContainer = $("#plots");
$.each(datasets, function(key, val) {
    choiceContainer.append(
		'<div class="plot_checkbox"><input type="checkbox" name="' + key + '" checked="checked" id="' + key + '"/><label for="'+ key + '">' + val.label + '</label></div>'
	);
});
choiceContainer.find("input").click(plotAccordingToChoices);

function addCommas(number) {
    number = '' + number;
    if (number.length > 3) {
        var mod = number.length % 3;
        var output = (mod > 0 ? (number.substring(0,mod)) : '');
        for (i=0 ; i < Math.floor(number.length / 3); i++) {
            if ((mod == 0) && (i == 0)) {
                output += number.substring(mod+ 3 * i, mod + 3 * i + 3);
			} else {
                output+= ',' + number.substring(mod + 3 * i, mod + 3 * i + 3);
			}
        }
        return (output);
    } else { return number; };
}

function showTooltip(x, y, contents) {
    $('<div id="tooltip">' + contents + '</div>').css( {
        position: 'absolute',
        display: 'none',
        top: y + 5,
        left: x + 5,
        border: '1px solid #fdd',
        padding: '2px',
        'background-color': '#fee',
        opacity: 0.80
    }).appendTo("body").fadeIn(200);
}

var previousPoint = null;
$("#placeholder").bind("plothover", function (event, pos, item) {
    $("#x").text(pos.x.toFixed(2));
    $("#y").text(pos.y.toFixed(2));

    if (item) {
        if (previousPoint != item.datapoint) {
            previousPoint = item.datapoint;

            $("#tooltip").remove();
            var x = item.datapoint[0].toFixed(2),
            y = item.datapoint[1].toFixed(2);

            showTooltip(item.pageX, item.pageY, item.series.label + " at " + x + " = " + y);
        }
    }
    else {
        $("#tooltip").remove();
        previousPoint = null;
    }
});

function plotAccordingToChoices() {
	var data = [];

	choiceContainer.find("input:checked").each(function () {
		var key = $(this).attr("name");
		if (key && datasets[key]) {
			data.push(datasets[key]);
		}
	});

	if (data.length > 0) {
		$.plot($("#placeholder"), data, {
				 legend: { show: true, position: "ne", margin: 20},
				 lines: { show: true },
				 points: { show: false },
				 xaxis: { ticks: 10 },
				 yaxis: { tickFormatter: function (v, axis) { return addCommas(v); }},
				 y2axis: { ticks: 10, min: 0, max: 100, tickFormatter: function (v, axis) { return v + "%";  } },
				 grid: { hoverable: true, clickable: true },
				 markings: [ { y2axis: { from:0, to: 100 } } ]
		});
	}
}
plotAccordingToChoices();
