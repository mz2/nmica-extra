package net.derkholm.nmica.extra.app;

import java.util.Arrays;
import java.util.List;

import org.bjv2.util.cli.CliTools;
import org.bjv2.util.cli.ConfigurationException;
import org.bjv2.util.cli.ConsoleMessages;
import org.bjv2.util.cli.UserLevel;

public class MotifAlignerApplicationDebug
{
  public static void main(String[] args)
          throws Throwable
  {
    MotifAligner app = new MotifAligner();
    List<String> argList = Arrays.asList(args);
    if (argList.indexOf("-help") >= 0) {
        ConsoleMessages.helpMessage(app, System.err, UserLevel.USER, 80);
        return;
    }
    try {
        args = CliTools.configureBean(app, args);
        app.main(args);
    } catch (ConfigurationException ex) {
        ConsoleMessages.errorMessage(app, System.err, 80, ex);
        System.exit(1);
    }
  }
}
