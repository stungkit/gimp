/* GIMP - The GNU Image Manipulation Program
 * Copyright (C) 1995 Spencer Kimball and Peter Mattis
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "config.h"
#include <gegl.h>

#include <mypaint-surface.h>

#include "paint-types.h"

#include "libgimpmath/gimpmath.h"

#include <cairo.h>
#include <gdk-pixbuf/gdk-pixbuf.h>
#include "libgimpcolor/gimpcolor.h"

#include "gimpmybrushoptions.h"
#include "gimpmybrushsurface.h"


struct _GimpMybrushSurface
{
  MyPaintSurface2     surface;
  GeglBuffer         *buffer;
  GeglBuffer         *paint_mask;
  gint                paint_mask_x;
  gint                paint_mask_y;
  gint                off_x;
  gint                off_y;
  GeglRectangle       dirty;
  GimpComponentMask   component_mask;
  GimpMybrushOptions *options;
};

/* --- Taken from mypaint-tiled-surface.c --- */
static inline float
calculate_rr (int   xp,
              int   yp,
              float x,
              float y,
              float aspect_ratio,
              float sn,
              float cs,
              float one_over_radius2)
{
    /* code duplication, see brush::count_dabs_to() */
    const float yy = (yp + 0.5f - y);
    const float xx = (xp + 0.5f - x);
    const float yyr=(yy*cs-xx*sn)*aspect_ratio;
    const float xxr=yy*sn+xx*cs;
    const float rr = (yyr*yyr + xxr*xxr) * one_over_radius2;
    /* rr is in range 0.0..1.0*sqrt(2) */
    return rr;
}

static inline float
calculate_r_sample (float x,
                    float y,
                    float aspect_ratio,
                    float sn,
                    float cs)
{
    const float yyr=(y*cs-x*sn)*aspect_ratio;
    const float xxr=y*sn+x*cs;
    const float r = (yyr*yyr + xxr*xxr);
    return r;
}

static inline float
sign_point_in_line (float px,
                    float py,
                    float vx,
                    float vy)
{
    return (px - vx) * (-vy) - (vx) * (py - vy);
}

static inline void
closest_point_to_line (float  lx,
                       float  ly,
                       float  px,
                       float  py,
                       float *ox,
                       float *oy)
{
    const float l2 = lx*lx + ly*ly;
    const float ltp_dot = px*lx + py*ly;
    const float t = ltp_dot / l2;
    *ox = lx * t;
    *oy = ly * t;
}


/* This works by taking the visibility at the nearest point
 * and dividing by 1.0 + delta.
 *
 * - nearest point: point where the dab has more influence
 * - farthest point: point at a fixed distance away from
 *                   the nearest point
 * - delta: how much occluded is the farthest point relative
 *          to the nearest point
 */
static inline float
calculate_rr_antialiased (int  xp,
                          int  yp,
                          float x,
                          float y,
                          float aspect_ratio,
                          float sn,
                          float cs,
                          float one_over_radius2,
                          float r_aa_start)
{
    /* calculate pixel position and borders in a way
     * that the dab's center is always at zero */
    float pixel_right = x - (float)xp;
    float pixel_bottom = y - (float)yp;
    float pixel_center_x = pixel_right - 0.5f;
    float pixel_center_y = pixel_bottom - 0.5f;
    float pixel_left = pixel_right - 1.0f;
    float pixel_top = pixel_bottom - 1.0f;

    float nearest_x, nearest_y; /* nearest to origin, but still inside pixel */
    float farthest_x, farthest_y; /* farthest from origin, but still inside pixel */
    float r_near, r_far, rr_near, rr_far;
    float center_sign, rad_area_1, visibilityNear, delta, delta2;

    /* Dab's center is inside pixel? */
    if( pixel_left<0 && pixel_right>0 &&
        pixel_top<0 && pixel_bottom>0 )
    {
        nearest_x = 0;
        nearest_y = 0;
        r_near = rr_near = 0;
    }
    else
    {
        closest_point_to_line( cs, sn, pixel_center_x, pixel_center_y, &nearest_x, &nearest_y );
        nearest_x = CLAMP( nearest_x, pixel_left, pixel_right );
        nearest_y = CLAMP( nearest_y, pixel_top, pixel_bottom );
        /* XXX: precision of "nearest" values could be improved
         * by intersecting the line that goes from nearest_x/Y to 0
         * with the pixel's borders here, however the improvements
         * would probably not justify the perdormance cost.
         */
        r_near = calculate_r_sample( nearest_x, nearest_y, aspect_ratio, sn, cs );
        rr_near = r_near * one_over_radius2;
    }

    /* out of dab's reach? */
    if( rr_near > 1.0f )
        return rr_near;

    /* check on which side of the dab's line is the pixel center */
    center_sign = sign_point_in_line( pixel_center_x, pixel_center_y, cs, -sn );

    /* radius of a circle with area=1
     *   A = pi * r * r
     *   r = sqrt(1/pi)
     */
    rad_area_1 = sqrtf( 1.0f / M_PI );

    /* center is below dab */
    if( center_sign < 0 )
    {
        farthest_x = nearest_x - sn*rad_area_1;
        farthest_y = nearest_y + cs*rad_area_1;
    }
    /* above dab */
    else
    {
        farthest_x = nearest_x + sn*rad_area_1;
        farthest_y = nearest_y - cs*rad_area_1;
    }

    r_far = calculate_r_sample( farthest_x, farthest_y, aspect_ratio, sn, cs );
    rr_far = r_far * one_over_radius2;

    /* check if we can skip heavier AA */
    if( r_far < r_aa_start )
        return (rr_far+rr_near) * 0.5f;

    /* calculate AA approximate */
    visibilityNear = 1.0f - rr_near;
    delta = rr_far - rr_near;
    delta2 = 1.0f + delta;
    visibilityNear /= delta2;

    return 1.0f - visibilityNear;
}
/* -- end mypaint code */

static inline float
calculate_alpha_for_rr (float rr,
                        float hardness,
                        float slope1,
                        float slope2)
{
  if (rr > 1.0f)
    return 0.0f;
  else if (rr <= hardness)
    return 1.0f + rr * slope1;
  else
    return rr * slope2 - slope2;
}

static GeglRectangle
calculate_dab_roi (float x,
                   float y,
                   float radius)
{
  int x0 = floor (x - radius);
  int x1 = ceil (x + radius);
  int y0 = floor (y - radius);
  int y1 = ceil (y + radius);

  return *GEGL_RECTANGLE (x0, y0, x1 - x0, y1 - y0);
}

static void
gimp_mypaint_surface_get_color_2 (MyPaintSurface2 *base_surface,
                                  gfloat           x,
                                  gfloat           y,
                                  gfloat           radius,
                                  gfloat          *color_r,
                                  gfloat          *color_g,
                                  gfloat          *color_b,
                                  gfloat          *color_a,
                                  gfloat           paint)
{
  GimpMybrushSurface *surface = (GimpMybrushSurface *)base_surface;
  GeglRectangle       dabRect;

  if (radius < 1.0f)
    radius = 1.0f;

  dabRect = calculate_dab_roi (x, y, radius);

  *color_r = 0.0f;
  *color_g = 0.0f;
  *color_b = 0.0f;
  *color_a = 0.0f;

  if (dabRect.width > 0 || dabRect.height > 0)
  {
    const float one_over_radius2 = 1.0f / (radius * radius);
    float sum_weight = 0.0f;
    float sum_r = 0.0f;
    float sum_g = 0.0f;
    float sum_b = 0.0f;
    float sum_a = 0.0f;

     /* Read in clamp mode to avoid transparency bleeding in at the edges */
    GeglBufferIterator *iter = gegl_buffer_iterator_new (surface->buffer, &dabRect, 0,
                                                         babl_format ("R'aG'aB'aA float"),
                                                         GEGL_BUFFER_READ,
                                                         GEGL_ABYSS_CLAMP, 2);
    if (surface->paint_mask)
      {
        GeglRectangle mask_roi = dabRect;
        mask_roi.x -= surface->paint_mask_x;
        mask_roi.y -= surface->paint_mask_y;
        gegl_buffer_iterator_add (iter, surface->paint_mask, &mask_roi, 0,
                                  babl_format ("Y float"),
                                  GEGL_ACCESS_READ, GEGL_ABYSS_NONE);
      }

    while (gegl_buffer_iterator_next (iter))
      {
        float *pixel = (float *)iter->items[0].data;
        float *mask;
        int iy, ix;

        if (surface->paint_mask)
          mask = iter->items[1].data;
        else
          mask = NULL;

        for (iy = iter->items[0].roi.y; iy < iter->items[0].roi.y + iter->items[0].roi.height; iy++)
          {
            float yy = (iy + 0.5f - y);
            for (ix = iter->items[0].roi.x; ix < iter->items[0].roi.x +  iter->items[0].roi.width; ix++)
              {
                /* pixel_weight == a standard dab with hardness = 0.5, aspect_ratio = 1.0, and angle = 0.0 */
                float xx = (ix + 0.5f - x);
                float rr = (yy * yy + xx * xx) * one_over_radius2;
                float pixel_weight = 0.0f;
                if (rr <= 1.0f)
                  pixel_weight = 1.0f - rr;
                if (mask)
                  pixel_weight *= *mask;

                sum_r += pixel_weight * pixel[RED];
                sum_g += pixel_weight * pixel[GREEN];
                sum_b += pixel_weight * pixel[BLUE];
                sum_a += pixel_weight * pixel[ALPHA];
                sum_weight += pixel_weight;

                pixel += 4;
                if (mask)
                  mask += 1;
              }
          }
      }

    if (sum_a > 0.0f && sum_weight > 0.0f)
      {
        sum_r /= sum_weight;
        sum_g /= sum_weight;
        sum_b /= sum_weight;
        sum_a /= sum_weight;

        sum_r /= sum_a;
        sum_g /= sum_a;
        sum_b /= sum_a;

        /* FIXME: Clamping is wrong because GEGL allows alpha > 1, this should probably re-multipy things */
        *color_r = CLAMP(sum_r, 0.0f, 1.0f);
        *color_g = CLAMP(sum_g, 0.0f, 1.0f);
        *color_b = CLAMP(sum_b, 0.0f, 1.0f);
        *color_a = CLAMP(sum_a, 0.0f, 1.0f);
      }
  }
}

static void
gimp_mypaint_surface_get_color_wrapper (MyPaintSurface *surface,
                                        float           x,
                                        float           y,
                                        float           radius,
                                        float          *color_r,
                                        float          *color_g,
                                        float          *color_b,
                                        float          *color_a)
{
  return gimp_mypaint_surface_get_color_2 ((MyPaintSurface2 *) surface, x, y, radius,
                                           color_r, color_g, color_b, color_a,
                                           -1.0);
}

static gint
gimp_mypaint_surface_draw_dab_2 (MyPaintSurface2 *base_surface,
                                 gfloat           x,
                                 gfloat           y,
                                 gfloat           radius,
                                 gfloat           color_r,
                                 gfloat           color_g,
                                 gfloat           color_b,
                                 gfloat           opaque,
                                 gfloat           hardness,
                                 gfloat           color_a,
                                 gfloat           aspect_ratio,
                                 gfloat           angle,
                                 gfloat           lock_alpha,
                                 gfloat           colorize,
                                 gfloat           posterize,
                                 gfloat           posterize_num,
                                 gfloat           paint)
{
  GimpMybrushSurface *surface = (GimpMybrushSurface *)base_surface;
  GeglBufferIterator *iter;
  GeglRectangle       dabRect;
  GimpComponentMask   component_mask  = surface->component_mask;
  /* XXX What spaces should we be working from and to? */
  const Babl         *rgb_to_hsl_fish = babl_fish (babl_format ("R'G'B' float"), babl_format ("HSL float"));
  const Babl         *hsl_to_rgb_fish = babl_fish (babl_format ("HSL float"), babl_format ("R'G'B' float"));

  const float one_over_radius2 = 1.0f / (radius * radius);
  const double angle_rad = angle / 360 * 2 * M_PI;
  const float cs = cos(angle_rad);
  const float sn = sin(angle_rad);
  float normal_mode;
  float segment1_slope;
  float segment2_slope;
  float r_aa_start;

  posterize     = CLAMP (posterize, 0.0f, 1.0f);
  posterize_num = CLAMP (ROUND (posterize_num * 100.0), 1, 128);
  paint         = CLAMP (paint, 0.0f, 1.0f);

  hardness = CLAMP (hardness, 0.0f, 1.0f);
  segment1_slope = -(1.0f / hardness - 1.0f);
  segment2_slope = -hardness / (1.0f - hardness);
  aspect_ratio = MAX (1.0f, aspect_ratio);

  r_aa_start = radius - 1.0f;
  r_aa_start = MAX (r_aa_start, 0);
  r_aa_start = (r_aa_start * r_aa_start) / aspect_ratio;

  normal_mode = opaque * (1.0f - colorize) * (1.0f - posterize);
  colorize = opaque * colorize;

  /* FIXME: This should use the real matrix values to trim aspect_ratio dabs */
  x += surface->off_x;
  y += surface->off_y;
  dabRect = calculate_dab_roi (x, y, radius);
  gegl_rectangle_intersect (&dabRect, &dabRect, gegl_buffer_get_extent (surface->buffer));

  if (dabRect.width <= 0 || dabRect.height <= 0)
    return 0;

  gegl_rectangle_bounding_box (&surface->dirty, &surface->dirty, &dabRect);

  iter = gegl_buffer_iterator_new (surface->buffer, &dabRect, 0,
                                   babl_format ("R'G'B'A float"),
                                   GEGL_BUFFER_READWRITE,
                                   GEGL_ABYSS_NONE, 2);
  if (surface->paint_mask)
    {
      GeglRectangle mask_roi = dabRect;
      mask_roi.x -= surface->paint_mask_x;
      mask_roi.y -= surface->paint_mask_y;
      gegl_buffer_iterator_add (iter, surface->paint_mask, &mask_roi, 0,
                                babl_format ("Y float"),
                                GEGL_ACCESS_READ, GEGL_ABYSS_NONE);
    }

  while (gegl_buffer_iterator_next (iter))
    {
      float *pixel = (float *)iter->items[0].data;
      float *mask;
      int iy, ix;

      if (surface->paint_mask)
        mask = iter->items[1].data;
      else
        mask = NULL;

      for (iy = iter->items[0].roi.y; iy < iter->items[0].roi.y + iter->items[0].roi.height; iy++)
        {
          for (ix = iter->items[0].roi.x; ix < iter->items[0].roi.x +  iter->items[0].roi.width; ix++)
            {
              float rr, base_alpha, alpha, dst_alpha, r, g, b, a;
              if (radius < 3.0f)
                rr = calculate_rr_antialiased (ix, iy, x, y, aspect_ratio, sn, cs, one_over_radius2, r_aa_start);
              else
                rr = calculate_rr (ix, iy, x, y, aspect_ratio, sn, cs, one_over_radius2);
              base_alpha = calculate_alpha_for_rr (rr, hardness, segment1_slope, segment2_slope);
              alpha = base_alpha * normal_mode;
              if (mask)
                alpha *= *mask;
              dst_alpha = pixel[ALPHA];
              /* a = alpha * color_a + dst_alpha * (1.0f - alpha);
               * which converts to: */
              a = alpha * (color_a - dst_alpha) + dst_alpha;
              r = pixel[RED];
              g = pixel[GREEN];
              b = pixel[BLUE];

              if (a > 0.0f)
                {
                  /* By definition the ratio between each color[] and pixel[] component in a non-pre-multipled blend always sums to 1.0f.
                   * Originally this would have been "(color[n] * alpha * color_a + pixel[n] * dst_alpha * (1.0f - alpha)) / a",
                   * instead we only calculate the cheaper term. */
                  float src_term = (alpha * color_a) / a;
                  float dst_term = 1.0f - src_term;
                  r = color_r * src_term + r * dst_term;
                  g = color_g * src_term + g * dst_term;
                  b = color_b * src_term + b * dst_term;
                }

              if (colorize > 0.0f && base_alpha > 0.0f)
                {
                  alpha = base_alpha * colorize;
                  a = alpha + dst_alpha - alpha * dst_alpha;
                  if (a > 0.0f)
                    {
                      float pixel_hsl[3], out_hsl[3];
                      float pixel_rgb[3] = {color_r, color_g, color_b};
                      float out_rgb[3]   = {r, g, b};
                      float src_term     = alpha / a;
                      float dst_term     = 1.0f - src_term;

                      /* Here I am completely unsure if the conversion are
                       * right, regarding color spaces. What is the color space
                       * of color_r/g/b arguments?
                       * TODO: this code should be double-checked.
                       */
                      babl_process (rgb_to_hsl_fish, pixel_rgb, pixel_hsl, 1);
                      babl_process (rgb_to_hsl_fish, out_rgb, out_hsl, 1);

                      out_hsl[0] = pixel_hsl[0];
                      out_hsl[1] = pixel_hsl[1];
                      babl_process (hsl_to_rgb_fish, out_hsl, out_rgb, 1);

                      r = (float)out_rgb[0] * src_term + r * dst_term;
                      g = (float)out_rgb[1] * src_term + g * dst_term;
                      b = (float)out_rgb[2] * src_term + b * dst_term;
                    }
                }

              if (posterize > 0.0f && base_alpha > 0.0f)
                {
                  alpha = base_alpha * posterize;
                  a     = alpha + dst_alpha - alpha * dst_alpha;
                  if (a > 0.0f)
                    {
                      gfloat post_pixel[3];
                      gfloat src_term = alpha / a;
                      gfloat dst_term = 1.0f - src_term;

                      post_pixel[0] = ROUND (r * posterize_num) / posterize_num;
                      post_pixel[1] = ROUND (g * posterize_num) / posterize_num;
                      post_pixel[2] = ROUND (b * posterize_num) / posterize_num;

                      r = post_pixel[0] * src_term + r * dst_term;
                      g = post_pixel[1] * src_term + g * dst_term;
                      b = post_pixel[2] * src_term + b * dst_term;
                    }
                }

              if (surface->options->no_erasing)
                a = MAX (a, pixel[ALPHA]);

              if (component_mask != GIMP_COMPONENT_MASK_ALL)
                {
                  if (component_mask & GIMP_COMPONENT_MASK_RED)
                    pixel[RED]   = r;
                  if (component_mask & GIMP_COMPONENT_MASK_GREEN)
                    pixel[GREEN] = g;
                  if (component_mask & GIMP_COMPONENT_MASK_BLUE)
                    pixel[BLUE]  = b;
                  if (component_mask & GIMP_COMPONENT_MASK_ALPHA)
                    pixel[ALPHA] = a;
                }
              else
                {
                  pixel[RED]   = r;
                  pixel[GREEN] = g;
                  pixel[BLUE]  = b;
                  pixel[ALPHA] = a;
                }

              pixel += 4;
              if (mask)
                mask += 1;
            }
        }
    }

  return 1;
}

static int
gimp_mypaint_surface_draw_dab_wrapper (MyPaintSurface *surface,
                                       float           x,
                                       float           y,
                                       float           radius,
                                       float           color_r,
                                       float           color_g,
                                       float           color_b,
                                       float           opaque,
                                       float           hardness,
                                       float           color_a,
                                       float           aspect_ratio,
                                       float           angle,
                                       float           lock_alpha,
                                       float           colorize)
{
  const gfloat posterize     = 0.0;
  const gfloat posterize_num = 1.0;
  const gfloat pigment       = 0.0;

  return gimp_mypaint_surface_draw_dab_2 ((MyPaintSurface2 *) surface, x, y, radius,
                                          color_r, color_g, color_b, opaque, hardness,
                                          color_a, aspect_ratio, angle, lock_alpha,
                                          colorize, posterize, posterize_num,
                                          pigment);
}

static void
gimp_mypaint_surface_begin_atomic (MyPaintSurface *base_surface)
{

}

static void
gimp_mypaint_surface_end_atomic_2 (MyPaintSurface2   *base_surface,
                                   MyPaintRectangles *rois)
{
  GimpMybrushSurface *surface = (GimpMybrushSurface *)base_surface;

  if (rois)
    {
      const gint roi_rects = rois->num_rectangles;

      for (gint i = 0; i < roi_rects; i++)
        {
          rois->rectangles[i].x         = surface->dirty.x;
          rois->rectangles[i].y         = surface->dirty.y;
          rois->rectangles[i].width     = surface->dirty.width;
          rois->rectangles[i].height    = surface->dirty.height;
          surface->dirty = *GEGL_RECTANGLE (0, 0, 0, 0);
        }
    }
}

static void
gimp_mypaint_surface_end_atomic_wrapper (MyPaintSurface   *surface,
                                         MyPaintRectangle *roi)
{
  MyPaintRectangles rois = {1, roi};

  gimp_mypaint_surface_end_atomic_2 ((MyPaintSurface2 *) surface, &rois);
}

static void
gimp_mypaint_surface_destroy (MyPaintSurface *base_surface)
{
  GimpMybrushSurface *surface = (GimpMybrushSurface *) base_surface;

  g_clear_object (&surface->buffer);
  g_clear_object (&surface->paint_mask);
  g_free (surface);
}

GimpMybrushSurface *
gimp_mypaint_surface_new (GeglBuffer         *buffer,
                          GimpComponentMask   component_mask,
                          GeglBuffer         *paint_mask,
                          gint                paint_mask_x,
                          gint                paint_mask_y,
                          GimpMybrushOptions *options)
{
  GimpMybrushSurface *surface = g_malloc0 (sizeof (GimpMybrushSurface));
  MyPaintSurface2    *s;

  mypaint_surface_init (&surface->surface.parent);
  s = &surface->surface;

  s->get_color_pigment    = gimp_mypaint_surface_get_color_2;
  s->draw_dab_pigment     = gimp_mypaint_surface_draw_dab_2;
  s->parent.begin_atomic  = gimp_mypaint_surface_begin_atomic;
  s->end_atomic_multi     = gimp_mypaint_surface_end_atomic_2;

  s->parent.draw_dab      = gimp_mypaint_surface_draw_dab_wrapper;
  s->parent.get_color     = gimp_mypaint_surface_get_color_wrapper;
  s->parent.end_atomic    = gimp_mypaint_surface_end_atomic_wrapper;

  s->parent.destroy       = gimp_mypaint_surface_destroy;

  surface->component_mask = component_mask;
  surface->options        = options;
  surface->buffer         = g_object_ref (buffer);
  if (paint_mask)
    surface->paint_mask   = g_object_ref (paint_mask);

  surface->paint_mask_x   = paint_mask_x;
  surface->paint_mask_y   = paint_mask_y;
  surface->dirty          = *GEGL_RECTANGLE (0, 0, 0, 0);

  surface->off_x          = 0;
  surface->off_y          = 0;

  return surface;
}

void
gimp_mypaint_surface_set_buffer (GimpMybrushSurface *surface,
                                 GeglBuffer         *buffer,
                                 gint                paint_mask_x,
                                 gint                paint_mask_y)
{
  g_object_unref (surface->buffer);

  surface->buffer = g_object_ref (buffer);

  surface->paint_mask_x = paint_mask_x;
  surface->paint_mask_y = paint_mask_y;
}

void
gimp_mypaint_surface_set_offset (GimpMybrushSurface *surface,
                                 gint                off_x,
                                 gint                off_y)
{
  surface->off_x = off_x;
  surface->off_y = off_y;
}

void
gimp_mypaint_surface_get_offset (GimpMybrushSurface *surface,
                                 gint               *off_x,
                                 gint               *off_y)
{
  *off_x = surface->off_x;
  *off_y = surface->off_y;
}
