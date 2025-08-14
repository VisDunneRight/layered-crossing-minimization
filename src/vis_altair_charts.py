import altair as alt


def draw_altair_scatter(data_points, x_axis, y_axis, color_field, x_title, y_title, chart_name, log_y_scale, plot_loess=False, loess_features=None, opacity=1, xdom=None):
    data = alt.Data(values=data_points)
    chart = alt.Chart(data).mark_circle(size=60).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title), scale=alt.Scale(domain=xdom) if xdom else alt.Scale()),
        y=alt.Y(f'{y_axis}:Q', scale=alt.Scale(type="log") if log_y_scale else alt.Scale(), axis=alt.Axis(title=y_title)) #,
        # color=alt.Color(f'{color_field}:N', scale=alt.Scale(scheme='dark2'))
    ).configure_mark(opacity=opacity)
    # .facet(column=f'{color_field}:N')
    # if plot_loess:
    #     for feature in loess_features:
    #         chart += chart.transform_filter(alt.FieldEqualPredicate(field=color_field, equal=feature)).transform_loess(
    #             x_axis, y_axis
    #         )
    chart.save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_line_chart(data_points, x_axis, y_axis, color_field, x_title, y_title, chart_name, log_y_scale):
    data = alt.Data(values=data_points)
    dom = list(set(dp[f"{color_field}"] for dp in data_points))
    # rng = ["#26547C", "#F0567A", "#E09D00", "#4ACB2A"]
    rng = ["#e15759", "#b07aa1", "#9c755f", "#f28e2b", "#ff9da7", "#4e79a7", "#59a14f", "#edc948", "#76b7b2"]
    chart = alt.Chart(data).mark_line(point={"filled": False, "fill": "white"}).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', scale=alt.Scale(type="log") if log_y_scale else None, axis=alt.Axis(title=y_title)),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(domain=dom, range=rng))
    )
    chart.save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_simple_line_chart(data_points, x_axis, y_axis, color_field, x_title, y_title, chart_name, xdom=None, ydom=None, coldom=None, rng2=False):
    data = alt.Data(values=data_points)
    dom = coldom if coldom is not None else sorted(list(set(dp[f"{color_field}"] for dp in data_points)))
    rng = ["#000000", "#e15759", "#4e79a7", "#b07aa1", "#f28e2b", "#9c755f", "#ff9da7", "#59a14f", "#4e79a7", "#ff9da7", "#59a14f"]  # "#edc948", "#76b7b2"]
    # color_palette = ["#ff8a80", "#ff1744", "#d50000", "#ea80fc", "#d500f9", "#aa00ff", "#82b1ff", "#2979ff", "#2962ff", "#ccff90", "#76ff03", "#64dd17", "#ffe57f", "#ffc400", "#ffab00"]
    color_palette = ["#e1bee7", "#ce93d8", "#ab47bc", "#7b1fa2", "#4a148c", "#b3e5fc", "#81d4fa", "#29b6f6", "#0288d1", "#01579b", "#c8e6c9", "#a5d6a7", "#66bb6a", "#388e3c", "#1b5e20", "#ffe0b2", "#ffcc80", "#ffa726", "#f57c00", "#e65100"]
    if rng2:
        rng = ["#000000", "#e15759", "#4e79a7", "#9c755f", "#ff9da7", "#59a14f"]
    chart = alt.Chart(data).mark_line(clip=True).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title), scale=alt.Scale(domain=xdom) if xdom else alt.Scale()),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title), scale=alt.Scale(domain=ydom) if ydom else alt.Scale()),
        # y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title, format='%'), scale=alt.Scale(domain=ydom) if ydom else alt.Scale()),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(domain=dom, range=color_palette)),
        # strokeDash=alt.StrokeDash("Dash:N", sort=["normal", "combined", "optimal"])
    )
    # .properties(
    #     width=800,
    #     height=600,
    #     autosize=alt.AutoSizeParams(
    #         type='fit',
    #         contains='padding'
    #     )
    # )
    chart.save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_scatter_with_regression_line(data_points, x_axis, y_axis, color_field, x_title, y_title, chart_name):
    data = alt.Data(values=data_points)
    chart = alt.Chart(data).mark_circle(size=60).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title)),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(scheme='dark2'))
    )
    chart2 = chart + chart.transform_regression(x_axis, y_axis).mark_line()
    chart2.save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_scatter_with_custom_line(scatter_data, line_data, x_axis, y_axis, color_field, x_title, y_title, chart_name):
    scatterdata = alt.Data(values=scatter_data)
    scatterchart = alt.Chart(scatterdata).mark_circle(size=60, opacity=0.9).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title)),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(scheme='dark2'))
    )
    linedata = alt.Data(values=line_data)
    linechart = alt.Chart(linedata).mark_line(color="black").encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title)),
        # color=alt.Color(f'{color_field}:N', scale=alt.Scale())
    )
    chart = scatterchart + linechart
    chart.save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_line_compare(data_points, x_axis, y_axis, facet_field, x_title, y_title, chart_name, log_y_scale, experiment_name):
    data = alt.Data(values=data_points)
    # dom = [f"{experiment_name}", "baseline"]
    rng = ["#26547C", "#F0567A"]
    chart = alt.Chart(data).mark_circle(size=60).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title)),
        # color=alt.Color(f'Technique:N', scale=alt.Scale(domain=dom, range=rng))
    )

    horizline = alt.Chart().mark_rule().encode(
        y='a:Q'
    )

    alt.layer(
        chart, horizline,
        data=data
    ).transform_calculate(
        a="100"
    ).facet(
        column=alt.Column(f'{facet_field}:N', sort=["junger_basic", "vertical_transitivity", "redundancy"])
    ).save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_colored_line_compare(data_points, x_axis, y_axis, facet_field, color_field, x_title, y_title, chart_name):
    data = alt.Data(values=data_points)
    # dom = [cat1name, cat2name]
    rng = ["#26547C", "#F0567A", "#E09D00"]
    chart = alt.Chart(data).mark_circle(size=60).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title)),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(range=rng))
    )

    horizline = alt.Chart().mark_rule().encode(
        y='a:Q'
    )

    alt.layer(
        chart, horizline,
        data=data
    ).transform_calculate(
        a="100"
    ).facet(
        column=alt.Column(f'{facet_field}:N', sort=["junger_basic", "vertical_transitivity"])
    ).save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})
