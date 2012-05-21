import org.openscience.cdk.Molecule;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.io.CMLReader;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import java.io.*;

public class Program
{
    public static void main(String[] args)
    {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileFilter(new javax.swing.filechooser.FileFilter()
        {
            public boolean accept(File f)
            {
                if (f.isDirectory()) return true;
                String s = f.getName();
                return s.length() > 4 && s.substring(s.length() - 4).equals(".cml");
            }
            
            public String getDescription()
            {
                return "Chemical Markup Language file (*.cml)";
            }
        });
        chooser.setAcceptAllFileFilterUsed(false);
        if (chooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
        {
            try
            {
                CMLReader reader = new CMLReader(new FileInputStream(chooser.getSelectedFile()));
                try
                {
                    ChemFile cf = new ChemFile();
                    cf = (ChemFile)reader.read(cf);
                    OrganicMolecule x = new OrganicMolecule(
                        new Molecule(ChemFileManipulator.getAllAtomContainers(cf).get(0)));
                    JOptionPane.showMessageDialog(null, x.getIUPACName(),
                        "Organic molecule name", JOptionPane.INFORMATION_MESSAGE);
                }
                catch (CDKException e2)
                {
                    JOptionPane.showMessageDialog(null, e2.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
                }
            }
            catch (FileNotFoundException e1)
            {
                JOptionPane.showMessageDialog(null, e1.getMessage(),
                    "Error", JOptionPane.ERROR_MESSAGE);
            }
        }
    }
}
