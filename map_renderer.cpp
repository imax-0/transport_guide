#include "map_renderer.h"
#include "sphere.h"
#include "sphere_projection.h"
#include "utils.h"
#include <cassert>
#include <algorithm>
#include <unordered_set>

using namespace std;

static Svg::Point ParsePoint(const Json::Node& json) {
  const auto& array = json.AsArray();
  return {
      array[0].AsDouble(),
      array[1].AsDouble()
  };
}

static Svg::Color ParseColor(const Json::Node& json) {
  if (json.IsString()) {
    return json.AsString();
  }
  const auto& array = json.AsArray();
  assert(array.size() == 3 || array.size() == 4);
  Svg::Rgb rgb{
      static_cast<uint8_t>(array[0].AsInt()),
      static_cast<uint8_t>(array[1].AsInt()),
      static_cast<uint8_t>(array[2].AsInt())
  };
  if (array.size() == 3) {
    return rgb;
  } else {
    return Svg::Rgba{rgb, array[3].AsDouble()};
  }
}

static vector<Svg::Color> ParseColors(const Json::Node& json) {
  const auto& array = json.AsArray();
  vector<Svg::Color> colors;
  colors.reserve(array.size());
  transform(begin(array), end(array), back_inserter(colors), ParseColor);
  return colors;
}

RenderSettings ParseRenderSettings(const Json::Dict& json) {
  RenderSettings result;
  result.max_width = json.at("width").AsDouble();
  result.max_height = json.at("height").AsDouble();
  result.padding = json.at("padding").AsDouble();
  result.outer_margin = json.at("outer_margin").AsDouble();
  result.palette = ParseColors(json.at("color_palette"));
  result.line_width = json.at("line_width").AsDouble();
  result.underlayer_color = ParseColor(json.at("underlayer_color"));
  result.underlayer_width = json.at("underlayer_width").AsDouble();
  result.stop_radius = json.at("stop_radius").AsDouble();
  result.bus_label_offset = ParsePoint(json.at("bus_label_offset"));
  result.bus_label_font_size = json.at("bus_label_font_size").AsInt();
  result.stop_label_offset = ParsePoint(json.at("stop_label_offset"));
  result.stop_label_font_size = json.at("stop_label_font_size").AsInt();

  const auto& layers_array = json.at("layers").AsArray();
  result.layers.reserve(layers_array.size());
  for (const auto& layer_node : layers_array) {
    result.layers.push_back(layer_node.AsString());
  }

  return result;
}

static unordered_set<string> FindReferenceStops(const Descriptions::BusesDict& buses_dict) {
	unordered_set<string> reference_stops;
	unordered_map<string, int> stops_rank;
	unordered_map<string, const Descriptions::Bus*> stop_bus;
	for(const auto& [_, bus_ptr] : buses_dict) {
		for(const auto& endpoint : bus_ptr->endpoints) {
			reference_stops.insert(endpoint);
		}

		for(const string& stop : bus_ptr->stops) {
			++stops_rank[stop];
			const auto [it, it_inserted] = stop_bus.emplace(stop, bus_ptr);
			if (!it_inserted && it->second != bus_ptr) {
				reference_stops.insert(stop);
			}
		}
	}

	for(auto& [name, rank] : stops_rank) {
		if (rank > 2) {
			reference_stops.insert(name);
		}
	}

	return reference_stops;
}

static unordered_map<string, Sphere::Point> RecalculateStopsDict(const Descriptions::StopsDict &stops_dict,
                                                                 const Descriptions::BusesDict& buses_dict) {
	unordered_set<string> reference_stops = FindReferenceStops(buses_dict);
	unordered_map<string, Sphere::Point> new_stops_dict;

	for(const auto& [temp, bus_ptr] : buses_dict) {
		size_t counter = 0;
		size_t i = 0;
		size_t j = 0;
		auto& stops = bus_ptr->stops;
		for (const string& stop : stops) {
			if (counter > 0) {
				auto it = reference_stops.find(stop);
				if (it != end(reference_stops)) {
					i = j;
					j = counter;

					double lon_step = (stops_dict.at(stops[j])->position.longitude - stops_dict.at(stops[i])->position.longitude) / (j - i);
					double lat_step = (stops_dict.at(stops[j])->position.latitude - stops_dict.at(stops[i])->position.latitude) / (j - i);

					for (size_t k = i; k <= j; ++k) {
						double new_lon_coord = stops_dict.at(stops[i])->position.longitude + lon_step * (k - i);
						double new_lat_coord = stops_dict.at(stops[i])->position.latitude + lat_step * (k - i);
						new_stops_dict[stops[k]] = { new_lat_coord, new_lon_coord };
					}
				}
			}
			++counter;
		}
	}

	for(const auto& [name, stop_ptr] : stops_dict) {
		if (new_stops_dict.find(name) == end(new_stops_dict)) {
			new_stops_dict[name] = stop_ptr->position;
		}
	}

	return new_stops_dict;
}

static map<string, Svg::Point> ComputeStopsCoords(const Descriptions::BusesDict& buses_dict,
		                                          const Descriptions::StopsDict& stops_dict,
                                                  const RenderSettings& render_settings) {

  unordered_map<string, Sphere::Point> new_stops_dict = RecalculateStopsDict(stops_dict, buses_dict);
  CoordsCompressor compressor(new_stops_dict, buses_dict);
  compressor.Compress(render_settings);

  map<string, Svg::Point> new_coords;
  for(const auto& [name, point] : new_stops_dict) {
	   new_coords[name] = { compressor.GetCoordOnMapLon(point.longitude),
			                compressor.GetCoordOnMapLat(point.latitude) };
  }

  return new_coords;
}

static unordered_map<string, Svg::Color> ChooseBusColors(const Descriptions::BusesDict& buses_dict,
                                                         const RenderSettings& render_settings) {
  const auto& palette = render_settings.palette;
  unordered_map<string, Svg::Color> bus_colors;
  int idx = 0;
  for (const auto& [bus_name, bus_ptr] : buses_dict) {
    bus_colors[bus_name] = palette[idx++ % palette.size()];
  }
  return bus_colors;
}

CoordsCompressor::CoordsCompressor(const unordered_map<string, Sphere::Point> &stops_dict,
                                   const Descriptions::BusesDict& buses_dict) {
	lons_.reserve(stops_dict.size());
	lats_.reserve(stops_dict.size());

	for (const auto& [_, position] : stops_dict) {
		lons_.push_back({position.longitude});
		lats_.push_back({position.latitude});
	}

	sort(begin(lons_), end(lons_));
	sort(begin(lats_), end(lats_));

	FillNeighbors(stops_dict, buses_dict);
}

void CoordsCompressor::FillNeighbors(const unordered_map<string, Sphere::Point>& stops_dict, const Descriptions::BusesDict& buses_dict) {

	for (const auto& [_, bus_ptr] : buses_dict) {
		auto& stops = bus_ptr->stops;
		for (size_t i = 1; i < stops.size(); ++i) {
			auto prev_stop = stops_dict.at(stops[i-1]);
			auto cur_stop = stops_dict.at(stops[i]);

			auto [min_lon, max_lon] = minmax(prev_stop.longitude, cur_stop.longitude);
			auto [min_lat, max_lat] = minmax(prev_stop.latitude, cur_stop.latitude);

			neighbour_lons_[max_lon].insert(min_lon);
			neighbour_lats_[max_lat].insert(min_lat);
		}
	}
}

int CoordsCompressor::ComputeMaxLonIdx() {
	for(size_t i = 1; i < lons_.size(); ++i) {
		auto& cur_lon = lons_[i];
		size_t j = i;
		int max_neighbour_idx = -1;
		while(j > 0) {
			--j;
			auto& prev_lon = lons_[j];
			if (neighbour_lons_[cur_lon.old_coord].find(prev_lon.old_coord) != neighbour_lons_[cur_lon.old_coord].end()) {
				if (prev_lon.idx > max_neighbour_idx) {
					max_neighbour_idx = prev_lon.idx;
				}
			}
		}
		cur_lon.idx = max_neighbour_idx + 1;
	}

	return max_element(begin(lons_), end(lons_), [](const auto& lhs, const auto& rhs) {
		return lhs.idx < rhs.idx;
	}) -> idx;
}

int CoordsCompressor::ComputeMaxLatIdx() {
	for(size_t i = 1; i < lats_.size(); ++i) {
		auto& cur_lat = lats_[i];
		size_t j = i;
		int max_neighbour_idx = -1;
		while(j > 0) {
			--j;
			auto& prev_lat = lats_[j];
			if (neighbour_lats_[cur_lat.old_coord].find(prev_lat.old_coord) != neighbour_lats_[cur_lat.old_coord].end()) {
				if (prev_lat.idx > max_neighbour_idx) {
					max_neighbour_idx = prev_lat.idx;
				}
			}
		}
		cur_lat.idx = max_neighbour_idx + 1;
	}
	return max_element(begin(lats_), end(lats_), [](const auto& lhs, const auto& rhs) {
		return lhs.idx < rhs.idx;
	}) -> idx;
}

void CoordsCompressor::Compress(const RenderSettings& render_settings) {
	if (lats_.empty() || lons_.empty()) {
		return;
	}
	double width = render_settings.max_width;
	double height = render_settings.max_height;
	double padding = render_settings.padding;

	size_t max_lon_idx = ComputeMaxLonIdx();
	double x_step = (max_lon_idx ? (width - 2 * padding) / max_lon_idx : 0);

	size_t max_lat_idx = ComputeMaxLatIdx();
	double y_step = (max_lat_idx ? (height - 2 * padding) / max_lat_idx : 0);

	{
		for (auto& lon : lons_) {
		  lon.coord_on_map = lon.idx * x_step + padding;
	  }
    }

	{
		for (auto& lat : lats_) {
		  lat.coord_on_map = height - padding - lat.idx * y_step;
	  }
    }
}

static map<string, Descriptions::Bus> CopyBusesDict(const Descriptions::BusesDict& source) {
  map<string, Descriptions::Bus> target;
  for (const auto& [name, data_ptr] : source) {
    target.emplace(name, *data_ptr);
  }
  return target;
}

MapRenderer::MapRenderer(const Descriptions::StopsDict& stops_dict,
                         const Descriptions::BusesDict& buses_dict,
                         const Json::Dict& render_settings_json)
    : render_settings_(ParseRenderSettings(render_settings_json)),
      buses_dict_(CopyBusesDict(buses_dict)),
      stops_coords_(ComputeStopsCoords(buses_dict, stops_dict, render_settings_)),
      bus_colors_(ChooseBusColors(buses_dict, render_settings_))
{
}

void MapRenderer::Print() {
	  for (auto& [bus_name, bus_ptr] : buses_dict_) {
		cout << bus_name << endl;
		for (const auto& stop : bus_ptr.stops) {
			cout << stop << " ";
		}
		cout << endl;
	  }
}

using WaitItem = TransportRouter::RouteInfo::WaitItem;
using BusItem = TransportRouter::RouteInfo::BusItem;

void MapRenderer::RenderMapBusLines(Svg::Document& svg) const {
  for (const auto& [bus_name, bus_ptr] : buses_dict_) {
    const auto& stops = bus_ptr.stops;
    if (stops.empty()) {
      continue;
    }
    Svg::Polyline line;
    line.SetStrokeColor(bus_colors_.at(bus_name))
        .SetStrokeWidth(render_settings_.line_width)
        .SetStrokeLineCap("round").SetStrokeLineJoin("round");
    for (const auto& stop_name : stops) {
      line.AddPoint(stops_coords_.at(stop_name));
    }
    svg.Add(line);
  }
}

void MapRenderer::RenderRouteBusLines(Svg::Document& svg, const TransportRouter::RouteInfo& route) const {
  for (auto& item : route.items) {
	if (holds_alternative<WaitItem>(item)) {
	  continue;
	}
	if (holds_alternative<BusItem>(item)) {
	  const auto& route_bus = get<BusItem>(item);
	  const string& bus_name = route_bus.bus_name;
	  const auto& stops = buses_dict_.at(bus_name).stops;

	  if (stops.empty()) {
	    continue;
	  }

	  Svg::Polyline line;
	  line.SetStrokeColor(bus_colors_.at(route_bus.bus_name))
	      .SetStrokeWidth(render_settings_.line_width)
	      .SetStrokeLineCap("round").SetStrokeLineJoin("round");

	  for (size_t i = route_bus.start_stop_idx_; i <= route_bus.finish_stop_idx_; ++i) {
		line.AddPoint(stops_coords_.at(stops[i]));
	  }
	  svg.Add(line);
	}
  }
}

void MapRenderer::RenderBusLabel(Svg::Document& svg, const string& bus_name, const string& stop_name) const {
  const auto& color = bus_colors_.at(bus_name);
  const auto point = stops_coords_.at(stop_name);

  const auto base_text =
      Svg::Text{}
      .SetPoint(point)
      .SetOffset(render_settings_.bus_label_offset)
      .SetFontSize(render_settings_.bus_label_font_size)
      .SetFontFamily("Verdana")
      .SetFontWeight("bold")
      .SetData(bus_name);
  svg.Add(
      Svg::Text(base_text)
      .SetFillColor(render_settings_.underlayer_color)
      .SetStrokeColor(render_settings_.underlayer_color)
      .SetStrokeWidth(render_settings_.underlayer_width)
      .SetStrokeLineCap("round").SetStrokeLineJoin("round")
  );
  svg.Add(
      Svg::Text(base_text)
      .SetFillColor(color)
  );
}

void MapRenderer::RenderMapBusLabels(Svg::Document& svg) const {
  for (const auto& [bus_name, bus_ptr] : buses_dict_) {
    const auto& stops = bus_ptr.stops;
    if (!stops.empty()) {
      for (const string& endpoint : bus_ptr.endpoints) {
        RenderBusLabel(svg, bus_name, endpoint);
      }
    }
  }
}

void MapRenderer::RenderRouteBusLabels(Svg::Document& svg, const TransportRouter::RouteInfo& route) const {
  for (auto& item : route.items) {
	if (holds_alternative<WaitItem>(item)) {
	  continue;
	}
	if (holds_alternative<BusItem>(item)) {
	  const auto& route_bus = get<BusItem>(item);
	  const string& bus_name = route_bus.bus_name;
	  const auto& stops = buses_dict_.at(bus_name).stops;

	  if (stops.empty()) {
	    continue;
	  }

	  auto& endpoints = buses_dict_.at(route_bus.bus_name).endpoints;

	  for (const size_t stop_idx : {route_bus.start_stop_idx_, route_bus.finish_stop_idx_}) {
		auto& stop = stops[stop_idx];
		if (stop_idx == 0
			|| stop_idx == stops.size() - 1
			|| find(endpoints.begin(), endpoints.end(), stop) != endpoints.end()) {
		  RenderBusLabel(svg, bus_name, stop);
		}
	  }
	}
  }
}

void MapRenderer::RenderStopPoint(Svg::Document& svg, Svg::Point stop_point) const {
  svg.Add(Svg::Circle{}
          .SetCenter(stop_point)
          .SetRadius(render_settings_.stop_radius)
          .SetFillColor("white"));
}

void MapRenderer::RenderMapStopPoints(Svg::Document& svg) const {
  for (const auto& [stop_name, stop_point] : stops_coords_) {
	RenderStopPoint(svg, stop_point);
  }
}

void MapRenderer::RenderRouteStopPoints(Svg::Document& svg, const TransportRouter::RouteInfo& route) const {
  for (auto& item : route.items) {
	if (holds_alternative<WaitItem>(item)) {
	  continue;
	}
	if (holds_alternative<BusItem>(item)) {
	  const auto& route_bus = get<BusItem>(item);
	  const string& bus_name = route_bus.bus_name;
	  const auto& stops = buses_dict_.at(bus_name).stops;

	  if (stops.empty()) {
	    continue;
	  }

	  for (size_t i = route_bus.start_stop_idx_; i <= route_bus.finish_stop_idx_; ++i) {
		RenderStopPoint(svg, stops_coords_.at(stops[i]));
	  }
	}
  }
}

void MapRenderer::RenderStopLabel(Svg::Document& svg, const string& stop_name, const Svg::Point stop_point) const {
  const auto base_text =
      Svg::Text{}
      .SetPoint(stop_point)
      .SetOffset(render_settings_.stop_label_offset)
      .SetFontSize(render_settings_.stop_label_font_size)
      .SetFontFamily("Verdana")
      .SetData(stop_name);
  svg.Add(
      Svg::Text(base_text)
      .SetFillColor(render_settings_.underlayer_color)
      .SetStrokeColor(render_settings_.underlayer_color)
      .SetStrokeWidth(render_settings_.underlayer_width)
      .SetStrokeLineCap("round").SetStrokeLineJoin("round")
  );
  svg.Add(
      Svg::Text(base_text)
      .SetFillColor("black")
  );
}

void MapRenderer::RenderMapStopLabels(Svg::Document& svg) const {
  for (const auto& [stop_name, stop_point] : stops_coords_) {
	RenderStopLabel(svg, stop_name, stop_point);
  }
}

void MapRenderer::RenderRouteStopLabels(Svg::Document& svg, const TransportRouter::RouteInfo& route) const {
  if (!route.items.empty()) {
    for (auto& item : route.items) {
	  if (holds_alternative<WaitItem>(item)) {
	    string wait_stop_name = get<WaitItem>(item).stop_name;
	    RenderStopLabel(svg, wait_stop_name, stops_coords_.at(wait_stop_name));
	  }
    }

    auto last_bus = get<BusItem>(route.items.back());
    const auto& stops = buses_dict_.at(last_bus.bus_name).stops;

    if (!stops.empty()) {
      const auto& finish_stop_name = stops[last_bus.finish_stop_idx_];
  	  RenderStopLabel(svg, finish_stop_name, stops_coords_.at(finish_stop_name));
    }
  }
}

const unordered_map<string, void (MapRenderer::*)(Svg::Document&) const> MapRenderer::MAP_LAYER_ACTIONS = {
    {"bus_lines",   &MapRenderer::RenderMapBusLines},
    {"bus_labels",  &MapRenderer::RenderMapBusLabels},
    {"stop_points", &MapRenderer::RenderMapStopPoints},
    {"stop_labels", &MapRenderer::RenderMapStopLabels},
};

const unordered_map<string, void (MapRenderer::*)(Svg::Document&, const TransportRouter::RouteInfo&) const> MapRenderer::ROUTE_LAYER_ACTIONS = {
    {"bus_lines",   &MapRenderer::RenderRouteBusLines},
    {"bus_labels",  &MapRenderer::RenderRouteBusLabels},
    {"stop_points", &MapRenderer::RenderRouteStopPoints},
    {"stop_labels", &MapRenderer::RenderRouteStopLabels},
};

Svg::Document MapRenderer::Render() const {
  Svg::Document svg;

  for (const auto& layer : render_settings_.layers) {
    (this->*MAP_LAYER_ACTIONS.at(layer))(svg);
  }

  return svg;
}

Svg::Document MapRenderer::RenderRoute(Svg::Document svg, const TransportRouter::RouteInfo& route) const {

  svg.Add(
    Svg::Rectangle {}
	.SetTopLeftPoint({-render_settings_.outer_margin, -render_settings_.outer_margin})
	.SetBottomRightPoint({render_settings_.outer_margin + render_settings_.max_width,
						  render_settings_.outer_margin + render_settings_.max_height})
	.SetFillColor(render_settings_.underlayer_color)
  );


  for (const auto& layer : render_settings_.layers) {
    (this->*ROUTE_LAYER_ACTIONS.at(layer))(svg, route);
  }

  return svg;
}
