classdef BBBDreaverApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        BBBDreaverUIFigure                matlab.ui.Figure
        StartDistLabel                  matlab.ui.control.Label
        startDistSpinner                matlab.ui.control.Spinner
        ExporttoExcelButton             matlab.ui.control.Button
        CreateanalysisobjectButton      matlab.ui.control.Button
        OpenREAVERGUIButton             matlab.ui.control.Button
        BrowseTestButton                matlab.ui.control.Button
        BrowseControlButton             matlab.ui.control.Button
        ProcessfoldersButton            matlab.ui.control.Button
        perivascularwidthpxSpinner      matlab.ui.control.Spinner
        perivascularwidthpxSpinnerLabel  matlab.ui.control.Label
        TestdirectoryEditField          matlab.ui.control.EditField
        TestdirectoryEditFieldLabel     matlab.ui.control.Label
        ControldirectoryEditField       matlab.ui.control.EditField
        ControldirectoryEditFieldLabel  matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            try
                evalin("base",'depndenciesImported')
            catch 
                DEV_INITIALIZE;
                assignin("base",'depndenciesImported',true)
            end
            
        end

        % Button pushed function: ProcessfoldersButton
        function ProcessfoldersButtonPushed(app, event)
            control_dir = app.ControldirectoryEditField.Value;
            test_dir = app.TestdirectoryEditField.Value;
            n_px = app.perivascularwidthpxSpinner.Value;
            diffu = app.startDistSpinner.Value;
            if ~isempty(control_dir)
                disp('Starting to work on control files')
                analyze_entire_folder(n_px,control_dir,diffu);
                disp('Finished processing control files')
            end
            if ~isempty(test_dir)
                disp('Starting to work on test files')
                analyze_entire_folder(n_px,test_dir,diffu); 
                disp('Finished processing test files')
            end
        end

        % Button pushed function: BrowseControlButton
        function BrowseControlButtonPushed(app, event)
            control_dir = uigetdir(pwd,'Choose control directory');
            app.ControldirectoryEditField.Value = control_dir;
            figure(app.BBBDreaverUIFigure)
        end

        % Button pushed function: BrowseTestButton
        function BrowseTestButtonPushed(app, event)
            test_dir = uigetdir(pwd,'Choose test directory');
            app.TestdirectoryEditField.Value = test_dir;
            figure(app.BBBDreaverUIFigure)
        end

        % Button pushed function: OpenREAVERGUIButton
        function OpenREAVERGUIButtonPushed(app, event)
            REAVER_GUI;
        end

        % Button pushed function: CreateanalysisobjectButton
        function CreateanalysisobjectButtonPushed(app, event)
            results = EB_analysis;
            assignin("base","results",results);
        end

        % Button pushed function: ExporttoExcelButton
        function ExporttoExcelButtonPushed(app, event)
            [control_filename,control_path] = ...
                uiputfile({'.csv','.xlsx'},...
                'Save control file','Control.csv');
            [test_filename,test_path] = ...
                uiputfile({'.csv','.xlsx'},...
                'Save test file','Test.csv');
            results = evalin('base','results');
            ths = 2:10;
            results.writecsv(ths,fullfile(control_path,control_filename),...
                fullfile(test_path,test_filename));
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create EBreaverUIFigure and hide until all components are created
            app.BBBDreaverUIFigure = uifigure('Visible', 'off');
            app.BBBDreaverUIFigure.Color = [0.9412 0.9412 0.9412];
            app.BBBDreaverUIFigure.Position = [100 100 535 368];
            app.BBBDreaverUIFigure.Name = 'EB reaver';

            % Create ControldirectoryEditFieldLabel
            app.ControldirectoryEditFieldLabel = uilabel(app.BBBDreaverUIFigure);
            app.ControldirectoryEditFieldLabel.HorizontalAlignment = 'right';
            app.ControldirectoryEditFieldLabel.Position = [70 254 94 22];
            app.ControldirectoryEditFieldLabel.Text = 'Control directory';

            % Create ControldirectoryEditField
            app.ControldirectoryEditField = uieditfield(app.BBBDreaverUIFigure, 'text');
            app.ControldirectoryEditField.Position = [179 253 204 23];

            % Create TestdirectoryEditFieldLabel
            app.TestdirectoryEditFieldLabel = uilabel(app.BBBDreaverUIFigure);
            app.TestdirectoryEditFieldLabel.HorizontalAlignment = 'right';
            app.TestdirectoryEditFieldLabel.Position = [70 212 77 22];
            app.TestdirectoryEditFieldLabel.Text = 'Test directory';

            % Create TestdirectoryEditField
            app.TestdirectoryEditField = uieditfield(app.BBBDreaverUIFigure, 'text');
            app.TestdirectoryEditField.Position = [179 211 204 23];

            % Create perivascularwidthpxSpinnerLabel
            app.perivascularwidthpxSpinnerLabel = uilabel(app.BBBDreaverUIFigure);
            app.perivascularwidthpxSpinnerLabel.HorizontalAlignment = 'right';
            app.perivascularwidthpxSpinnerLabel.Position = [70 160 124 22];
            app.perivascularwidthpxSpinnerLabel.Text = 'perivascular width [px]';

            % Create perivascularwidthpxSpinner
            app.perivascularwidthpxSpinner = uispinner(app.BBBDreaverUIFigure);
            app.perivascularwidthpxSpinner.Position = [209 155 56 31];

            % Create ProcessfoldersButton
            app.ProcessfoldersButton = uibutton(app.BBBDreaverUIFigure, 'push');
            app.ProcessfoldersButton.ButtonPushedFcn = createCallbackFcn(app, @ProcessfoldersButtonPushed, true);
            app.ProcessfoldersButton.Position = [283 137 122 30];
            app.ProcessfoldersButton.Text = 'Process folders';

            % Create BrowseControlButton
            app.BrowseControlButton = uibutton(app.BBBDreaverUIFigure, 'push');
            app.BrowseControlButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseControlButtonPushed, true);
            app.BrowseControlButton.Position = [404 253 72 22];
            app.BrowseControlButton.Text = 'Browse';

            % Create BrowseTestButton
            app.BrowseTestButton = uibutton(app.BBBDreaverUIFigure, 'push');
            app.BrowseTestButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseTestButtonPushed, true);
            app.BrowseTestButton.Position = [404 212 72 22];
            app.BrowseTestButton.Text = 'Browse';

            % Create OpenREAVERGUIButton
            app.OpenREAVERGUIButton = uibutton(app.BBBDreaverUIFigure, 'push');
            app.OpenREAVERGUIButton.ButtonPushedFcn = createCallbackFcn(app, @OpenREAVERGUIButtonPushed, true);
            app.OpenREAVERGUIButton.Position = [173 306 215 34];
            app.OpenREAVERGUIButton.Text = 'Open REAVER GUI for segmentation';

            % Create CreateanalysisobjectButton
            app.CreateanalysisobjectButton = uibutton(app.BBBDreaverUIFigure, 'push');
            app.CreateanalysisobjectButton.ButtonPushedFcn = createCallbackFcn(app, @CreateanalysisobjectButtonPushed, true);
            app.CreateanalysisobjectButton.Position = [181 60 208 25];
            app.CreateanalysisobjectButton.Text = 'Create analysis object';


            % Create ExporttoExcelButton
            app.ExporttoExcelButton = uibutton(app.BBBDreaverUIFigure, 'push');
            app.ExporttoExcelButton.ButtonPushedFcn = createCallbackFcn(app, @ExporttoExcelButtonPushed, true);
            app.ExporttoExcelButton.Position = [221 21 121 27];
            app.ExporttoExcelButton.Text = 'Export to Excel';

            % Create startDistSpinner
            app.startDistSpinner = uispinner(app.BBBDreaverUIFigure);
            app.startDistSpinner.Position = [209 117 57 31];

            % Create StartDistLabel
            app.StartDistLabel = uilabel(app.BBBDreaverUIFigure);
            app.StartDistLabel.HorizontalAlignment = 'right';
            app.StartDistLabel.Position = [31 121 164 22];
            app.StartDistLabel.Text = 'Distance from vessel wall [px]';

            % Show the figure after all components are created
            app.BBBDreaverUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = EBreaverApp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.EBreaverUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.BBBDreaverUIFigure)
        end
    end
end