package org.panda.gem;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.svggen.SVGGraphics2DIOException;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * Creates an SVG figure showing GEM results.
 *
 * @author Emek Demir
 * @author Ozgun Babur
 */
public class ModPrint
{
	/**
	 * Flag to control grouping method. Results can be grouped by modulator (groupMod = true) or by
	 * targets (groupMod = false).
	 */
	static boolean groupMod = true;

	public static final AlphaComposite text =
		AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 1f);
	public static final AlphaComposite box =
		AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f);
	public static final Font modfont = new Font("Verdana", Font.PLAIN, 12);
	public static final Font targetfont = new Font("Verdana", Font.PLAIN, 8);
	public static final Font labelfont = new Font("Verdana", Font.PLAIN, 6);

	/**
	 * Each figure is drawn for just one factor gene.
	 */
	protected String factorName;
	private HashMap<String, Group> model;

	//public static Font upFont;
	//public static Font downFont;

	public static void main(String[] args) throws Throwable
	{
		String file = "result/Result_fdr0.05_var1.0_AR_expo_prostate";
		new ModPrint().generateGEMPlot(file + ".xls", file + ".svg");
	}

	/**
	 * Draws the list of group.
	 *
	 * @param g2d
	 */
	public void paint(Graphics2D g2d)
	{
		List<Group> modlist= new ArrayList<Group>(model.values());
		Collections.sort(modlist);

		int i = 0;
		Group[][] page = new Group[(modlist.size() / 5) + (modlist.size() % 5 == 0 ? 0 : 1)][5];
		for (Group mod : modlist)
		{
			page[i / 5][i % 5] = mod;
			i++;
		}
		int rowStart = 60;
		for (int rowi = 0; rowi < page.length; rowi++)
		{
			Group[] row = page[rowi];
			int top = 0;
			int bottom = 0;
			for (Group mod : row)
			{
				if (mod != null)
				{
					top = Math.max(top, mod.getTop());
					bottom = Math.max(bottom, mod.getBottom());
				}
			}

			int rowCenter = rowStart + 10 * top + 6;
			int rowbottom = rowCenter + 10 * bottom + 42;
			drawFactorLine(g2d, rowCenter,0);

			for (int col = 0; col < row.length; col++)
			{
				Group mod = row[col];
				if (mod != null)
				{
					drawGroup(g2d, mod, 50 + 120 * col, rowStart, rowCenter, rowbottom);
				}
			}
			rowStart = rowbottom;
		}
	}

	protected void drawFactorLine(Graphics2D g2d, int rowCenter, int offset)
	{
		setText(g2d);
		g2d.setFont(modfont);
		g2d.drawString(factorName, offset+10, rowCenter + 5);
		g2d.setFont(labelfont);
		g2d.drawString("activates", offset+10, rowCenter - 7);
		g2d.setFont(labelfont);
		g2d.drawString("inhibits", offset+10, rowCenter + 12);
	}

	/**
	 * Draws a group.
	 *
	 * @param g2d
	 * @param mod
	 * @param x
	 * @param rowStart
	 * @param rowCenter
	 * @param rowBottom
	 */
	protected void drawGroup(Graphics2D g2d, Group mod, int x, int rowStart, int rowCenter, int rowBottom)
	{
		mod.sortMembers();

		g2d.setColor(Color.LIGHT_GRAY);
		g2d.drawRect(x + 27, rowStart - 25, 110, rowBottom - rowStart - 10);

		setText(g2d);
		g2d.setFont(modfont);
		String name = mod.getName();
		FontMetrics fm = g2d.getFontMetrics();
		int width = fm.stringWidth(name);
		int margin = (120 - width) / 2 + 27;

		g2d.drawString(name, x + margin, rowStart - 10);
		g2d.setFont(labelfont);
		g2d.drawString("Enhances", x + 33, rowStart);
		g2d.drawString("Attenuates", x + 67, rowStart);
		g2d.drawString("Inverts", x + 107, rowStart);
		for (Cat cat : Cat.values())
		{
			List<String> targets = mod.get(cat);
			if (targets != null)
			{
				setBox(g2d, cat);
				Rectangle bounds = cat.getBounds(mod.getHeight(cat) * 10, x, rowCenter);
				g2d.fill(bounds);
				setText(g2d);
				int i = 0;
				for (String target : targets)
				{
					i++;
					drawTarget(g2d, bounds.x + 2, bounds.y - 2 + i * 10, target);
				}
			}
		}
	}

	protected void drawTarget(Graphics2D g2d, int x, int y, String target)
	{
		g2d.drawString(target, x, y);
	}

	/**
	 * Sets color to the category.
	 *
	 * @param g2d
	 * @param cat
	 */
	private void setBox(Graphics2D g2d, Cat cat)
	{
		g2d.setPaint(cat.getColor());
		g2d.setComposite(box);
	}

	/**
	 * Switches to text color.
	 *
	 * @param g2d
	 */
	private void setText(Graphics2D g2d)
	{
		g2d.setPaint(Color.BLACK);
		g2d.setComposite(text);
	}


	public void generateGEMPlot(String tripFile, String outFile) throws IOException
	{
		generateGEMPlot(Triplet.load(tripFile), outFile);
	}

	public void generateGEMPlot(List<Triplet> trips, String outFile) throws IOException
	{
		// Prepare model
		process(trips);

		// Create an instance of the SVG Generator.
		generateSVG(outFile);
	}

	/**
	 * Render using SVG Graphics2D implementation.
	 *
	 * @param outFile
	 * @throws UnsupportedEncodingException
	 * @throws FileNotFoundException
	 * @throws SVGGraphics2DIOException
	 */
	private void generateSVG(String outFile)
		throws UnsupportedEncodingException, FileNotFoundException, SVGGraphics2DIOException
	{
		// Get a DOMImplementation.
		DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();

		// Create an instance of org.w3c.dom.Document.
		String svgNS = "http://www.w3.org/2000/svg";
		Document document = domImpl.createDocument(svgNS, "svg", null);
		SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
		paint(svgGenerator);

		// Finally, stream out SVG to the standard output using
		// UTF-8 encoding.
		boolean useCSS = true; // we want to use CSS style attributes
		Writer out = new OutputStreamWriter(new FileOutputStream(outFile), "UTF-8");
		svgGenerator.stream(out, useCSS);
	}

	/**
	 * Reads triplets and generates Mod list.
	 *
	 * @param trips
	 * @return
	 * @throws IOException
	 */
	protected void process(List<Triplet> trips) throws IOException
	{
		this.model = new HashMap<>();

		for (Triplet t : trips)
		{
			if (factorName == null)
			{
				factorName = t.F.symbol;
			}

			String modul = groupMod ? t.M.symbol : t.T.symbol;

			if (!model.containsKey(modul))
			{
				model.put(modul, new Group(modul));
			}
			Group mod = model.get(modul);

			String tar = groupMod ? t.T.symbol : t.M.symbol;

			if (t.cat == null) continue;

			Cat cat = Cat.valueOf(t.cat.toString());

			if (mod.get(cat) == null)
			{
				mod.put(cat, new ArrayList<>());
			}
			List<String> tars = mod.get(cat);

			if (!tars.contains(tar))
			{
				tars.add(tar);
			}
		}
//		filterGroups();
	}

	/**
	 * Filters out unqualifying groups.
	 *	 */
	private void filterGroups()
	{
		for (String key : new HashSet<>(model.keySet()))
		{
			Group mod = model.get(key);

			if (mod.getSize() < 3)
			{
				model.remove(key);
			}
		}
	}

	/**
	 * This can be a group of targets with a common modulator, or a group of modulators with a common
	 * target.
	 */
	static class Group implements Comparable
	{
		private String name;
		HashMap<Cat, List<String>> cats = new HashMap<>();

		public Group(String modName)
		{
			this.name = modName;
		}

		public int compareTo(Object o)
		{
			return ((Group) o).getSize() - this.getSize();
		}

		public int getSize()
		{
			int i = 0;
			for (List<String> strings : cats.values())
			{
				i += strings == null ? 0 : strings.size();
			}

			return i;
		}

		public int getTop()
		{
			return Math.max(getHeight(Cat.ENHANCES_ACTIVATION),
				Math.max(getHeight(Cat.ATTENUATES_ACTIVATION),
					getHeight(Cat.INVERTS_ACTIVATION)));
		}

		public int getBottom()
		{
			return Math.max(getHeight(Cat.ENHANCES_INHIBITION),
				Math.max(getHeight(Cat.ATTENUATES_INHIBITION),
					getHeight(Cat.INVERTS_INHIBITION)));
		}

		public int getHeight(Cat cat)
		{
			return cats.get(cat) == null ? 0 : cats.get(cat).size();
		}

		public String getName()
		{
			return name;
		}

		public List<String> get(Cat cat)
		{
			return cats.get(cat);
		}

		public void put(Cat cat, List<String> targets)
		{
			this.cats.put(cat, targets);
		}

		public void sortMembers()
		{
			cats.values().forEach(Collections::sort);
		}
	}

	/**
	 * A categoty of modulation.
	 */
	enum Cat
	{
		ENHANCES_ACTIVATION(33, 137, 33, 1, 1),
		ATTENUATES_ACTIVATION(144, 238, 144, 1, 2),
		INVERTS_ACTIVATION(253, 48, 48, 1, 3),
		ENHANCES_INHIBITION(106, 33, 137, 0, 1),
		ATTENUATES_INHIBITION(210, 143, 238, 0, 2),
		INVERTS_INHIBITION(253, 207, 47, 0, 3);

		private Color color;
		private int v, h;

		Cat(int r, int g, int b, int v, int h)
		{
			color = new Color(r, g, b);
			this.v = v;
			this.h = h;
		}

		public Color getColor()
		{
			return color;
		}

		public Rectangle getBounds(int height, int x, int y)
		{
			return new Rectangle(x + h * 34, y - v * (height + 4) + 2, 32, height);
		}
	}
}
