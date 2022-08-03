#pragma once

#include "descriptions.h"
#include "json.h"
#include "svg.h"
#include "transport_router.h"

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>


struct RenderSettings {
  double max_width;
  double max_height;
  double padding;
  double outer_margin;
  std::vector<Svg::Color> palette;
  double line_width;
  Svg::Color underlayer_color;
  double underlayer_width;
  double stop_radius;
  Svg::Point bus_label_offset;
  int bus_label_font_size;
  Svg::Point stop_label_offset;
  int stop_label_font_size;
  std::vector<std::string> layers;
};

class CoordsCompressor {
public:
	using ReferenceStopsIdx = int;

	CoordsCompressor(const std::unordered_map<std::string, Sphere::Point>& stops_dict,
			         const Descriptions::BusesDict& buses_dict);
	void FillNeighbors(const std::unordered_map<std::string, Sphere::Point>& stops_dict, const Descriptions::BusesDict& buses_dict);
	void Compress(const RenderSettings& render_settings);

	double GetCoordOnMapLon(double old_coord) {
		return lower_bound(begin(lons_), end(lons_), CoordsInfo{old_coord})->coord_on_map;
	}
	double GetCoordOnMapLat(double old_coord) {
		return lower_bound(begin(lats_), end(lats_), CoordsInfo{old_coord})->coord_on_map;
	}

private:
	std::unordered_map<double, std::unordered_set<double>> neighbour_lons_;
	std::unordered_map<double, std::unordered_set<double>> neighbour_lats_;

	struct CoordsInfo {
		double old_coord;
		double coord_on_map = 0;
		int idx = 0;

		bool operator< (const CoordsInfo& coordinate) const {
			return old_coord < coordinate.old_coord;
		}
	};

	std::vector<CoordsInfo> lons_;
	std::vector<CoordsInfo> lats_;

	int ComputeMaxLonIdx();
	int ComputeMaxLatIdx();
};

class MapRenderer {
public:
  MapRenderer(const Descriptions::StopsDict& stops_dict,
              const Descriptions::BusesDict& buses_dict,
              const Json::Dict& render_settings_json);

  Svg::Document Render() const;
  Svg::Document RenderRoute(Svg::Document svg, const TransportRouter::RouteInfo& route) const;
  void Print();

private:
  RenderSettings render_settings_;
  std::map<std::string, Descriptions::Bus> buses_dict_;
  std::map<std::string, Svg::Point> stops_coords_;
  std::unordered_map<std::string, Svg::Color> bus_colors_;

  void RenderBusLabel(Svg::Document& svg, const std::string& bus_name, const std::string& stop_name) const;
  void RenderStopPoint(Svg::Document& svg, const Svg::Point stop_point) const;
  void RenderStopLabel(Svg::Document& svg, const std::string& stop_name, const Svg::Point stop_point) const;

  void RenderMapBusLines(Svg::Document& svg) const;
  void RenderMapBusLabels(Svg::Document& svg) const;
  void RenderMapStopPoints(Svg::Document& svg) const;
  void RenderMapStopLabels(Svg::Document& svg) const;

  void RenderRouteBusLines(Svg::Document& svg, const TransportRouter::RouteInfo& route) const;
  void RenderRouteBusLabels(Svg::Document& svg, const TransportRouter::RouteInfo& route) const;
  void RenderRouteStopPoints(Svg::Document& svg, const TransportRouter::RouteInfo& route) const;
  void RenderRouteStopLabels(Svg::Document& svg, const TransportRouter::RouteInfo& route) const;

  static const std::unordered_map<std::string, void (MapRenderer::*)(Svg::Document&) const> MAP_LAYER_ACTIONS;
  static const std::unordered_map<std::string, void (MapRenderer::*)(Svg::Document&, const TransportRouter::RouteInfo&) const> ROUTE_LAYER_ACTIONS;
};
