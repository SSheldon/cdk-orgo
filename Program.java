import org.openscience.cdk.Molecule;
import org.openscience.cdk.templates.MoleculeFactory;
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
                    JOptionPane.showMessageDialog(null, x.getIUPACName());
                }
                catch (CDKException e2)
                {
                    System.out.println(e2);
                }
            }
            catch (FileNotFoundException e1)
            {
                System.out.println(e1);
            }
        }
    }
}
