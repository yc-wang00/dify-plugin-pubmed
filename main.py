from dify_plugin import Plugin, DifyPluginEnv

# Set a reasonable timeout for potentially slow PubMed requests
plugin = Plugin(DifyPluginEnv(MAX_REQUEST_TIMEOUT=120))

if __name__ == '__main__':
    plugin.run()
